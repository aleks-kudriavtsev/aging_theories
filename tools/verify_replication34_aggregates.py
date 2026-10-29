#!/usr/bin/env python3
"""Independent stdlib mathematics for replication34 audit (no producer imports).

Prepared from frozen4df636ea protocol before production outputs. This module
does not simulate cases, fit models, read participants or make network calls.
Final exported-result/schema adapter is added after the runner is available.
"""
from fractions import Fraction
from functools import lru_cache
import argparse
from collections import Counter
import hashlib
import json
import math
from pathlib import Path
import re
from statistics import NormalDist
import subprocess
from urllib.parse import urlparse

PROTOCOL_COMMIT='4df636ea18ca8bb29b05788c37800438566426ec'

def close(a,b,rel=1e-10,absolute=1e-12):
    return (a is None and b is None) or (type(a) in (int,float) and type(b) in (int,float)
        and math.isfinite(a) and math.isfinite(b) and math.isclose(a,b,rel_tol=rel,abs_tol=absolute))

@lru_cache(maxsize=32)
def binomial_probabilities(n,p):
    """Independent normalized relative-PMF recurrence, no scipy/binomial helper.

    Frozen scenarios use n<=20000; the resource bound is an auditor budget,
    not a proposed clinical sample-size restriction.
    """
    if type(n) is not int or not 0<=n<=1000000 or not 0<=p<=1:
        raise ValueError('Binomial audit input outside the declared calculation budget')
    masses=[0.]*(n+1)
    if p==0:masses[0]=1.;return tuple(masses)
    if p==1:masses[n]=1.;return tuple(masses)
    mode=min(n,math.floor((n+1)*p));masses[mode]=1.
    for k in range(mode,0,-1):masses[k-1]=masses[k]*k/(n-k+1)*(1-p)/p
    for k in range(mode,n):masses[k+1]=masses[k]*(n-k)/(k+1)*p/(1-p)
    total=math.fsum(masses)
    return tuple(value/total for value in masses)

def binomial_cdf(k,n,p):
    if k<0:return 0.
    if k>=n:return 1.
    return math.fsum(binomial_probabilities(n,p)[:k+1])

def binomial_quantile(q,n,p):
    if not 0<q<1:raise ValueError('Audit inverse CDF uses open quantile levels')
    if p==0:return 0
    if p==1:return n
    total=correction=0.
    for k,mass in enumerate(binomial_probabilities(n,p)):
        increment=mass-correction;new=total+increment
        correction=(new-total)-increment;total=new
        if total>=q:return k
    return n

def expected_plan(s):
    """Alternative algebra from explicit locked inputs, never33B estimates."""
    mu=s['mu_d'];bound=max(0.,1-mu*mu)
    if 'sd_d' in s:
        variance=s['sd_d']**2
        variance_source='specified_individual_paired_loss_SD'
    else:
        variance=(s['se_reference']**2)*(s['n_reference']/s['deff_reference'])
        variance_source='anchored_SE_inverted_under_conditional_loss_variance_assumption'
    clipped=variance>bound;variance=min(variance,bound)
    z=NormalDist().inv_cdf(1-s['alpha']/2)
    se=math.sqrt(variance*s['deff_d']/s['n_evaluable'])
    n_precision=max(1,math.ceil((z*z*s['deff_d']*variance)/(s['half_width']**2)))
    availability=s['q_joint']*s['q_known_conditional']
    n_recruit=max(1,math.ceil(n_precision/availability)) if availability>0 else None
    n=s['n_evaluable'];p=s['raw_evaluable_event_probability'];floor=s.get('technical_event_floor')
    independent=s['event_count_model']=='independent_binomial'
    return dict(mean_d=mu,paired_loss_variance=variance,paired_loss_SD=math.sqrt(variance),
        variance_source=variance_source,moment_variance_roundoff_clipped=clipped,
        normal_quantile=z,conditional_mean_variance=variance*s['deff_d']/n,
        conditional_SE=se,conditional_normal_half_width=z*se,
        conditional_n_for_half_width=n_precision,joint_evaluable_yield_probability=availability,
        candidate_recruitment_expectation=n_recruit,degenerate_zero_loss_variance=variance==0,
        event_expected=n*p,event_variance=s['deff_event']*n*p*(1-p),
        event_quantiles=[binomial_quantile(s['alpha']/2,n,p),binomial_quantile(1-s['alpha']/2,n,p)] if independent else None,
        event_probability_below_floor=binomial_cdf(floor-1,n,p) if independent and floor is not None else None)

def prediction_atom_moments(atoms):
    """Exact decimal-rational sums using literal squared losses, not factoring."""
    f=lambda value:Fraction(str(value))
    qsum=sum(f(a['q']) for a in atoms)
    rows=[]
    for atom in atoms:
        q=f(atom['q'])/qsum;r=f(atom['r']);new=f(atom['p_new']);ref=f(atom['p_ref'])
        d0=new**2-ref**2;d1=(1-new)**2-(1-ref)**2
        rows.append((q,r,d0,d1))
    pi=sum(q*r for q,r,d0,d1 in rows)
    mu=sum(q*((1-r)*d0+r*d1) for q,r,d0,d1 in rows)
    second=sum(q*((1-r)*d0**2+r*d1**2) for q,r,d0,d1 in rows)
    variance=second-mu**2
    within=sum(q*(((1-r)*d0**2+r*d1**2)-((1-r)*d0+r*d1)**2) for q,r,d0,d1 in rows)
    between=sum(q*(((1-r)*d0+r*d1)-mu)**2 for q,r,d0,d1 in rows)
    fourth=sum(q*((1-r)*(d0-mu)**4+r*(d1-mu)**4) for q,r,d0,d1 in rows)
    return dict(input_probability_sum=float(qsum),raw_event_probability=float(pi),mean_d=float(mu),
        second_moment_d=float(second),variance_d=float(variance),SD_d=math.sqrt(float(variance)),
        mean_conditional_Y_variance=float(within),variance_of_conditional_means=float(between),
        fourth_central_moment_d=float(fourth))

def mean_fourth_moment(individual_variance,individual_fourth,j,m,rho):
    """Exact shared-pair-mixture cluster moment, then mean of j clusters."""
    v=individual_variance
    cluster_variance=v*(rho+(1-rho)/m)
    cluster_fourth=rho*individual_fourth+(1-rho)*(individual_fourth+3*(m-1)*v*v)/(m**3)
    mean_variance=cluster_variance/j
    mean_fourth=(cluster_fourth+3*(j-1)*cluster_variance**2)/(j**3)
    return mean_variance,mean_fourth

def expected_simulation(case):
    exact=prediction_atom_moments(case['atoms']);j=case['n_full_clusters'];m=case['cluster_size'];rho=case['rho'];n=j*m
    deff=1+(m-1)*rho;v,mu4=mean_fourth_moment(exact['variance_d'],exact['fourth_central_moment_d'],j,m,rho)
    z=NormalDist().inv_cdf(1-case['alpha']/2);p=exact['raw_event_probability'];r=case['replicates']
    vy=p*(1-p);fourth_y=vy*((1-p)**3+p**3)
    mean_y_v,mean_y_fourth=mean_fourth_moment(vy,fourth_y,j,m,rho)
    count_variance=(n*n)*mean_y_v;count_fourth=(n**4)*mean_y_fourth
    return dict(exact_moments=exact,n_participants=n,deff=deff,exact_mean_d=exact['mean_d'],
        exact_variance_mean_d=v,exact_normal_half_width=z*math.sqrt(v),
        MCSE_mean_d=math.sqrt(v/r),exact_expected_event_count=n*p,
        exact_event_count_variance_in_this_toy=count_variance,MCSE_mean_event_count=math.sqrt(count_variance/r),
        theoretical_fourth_central_moment_mean_d=mu4,
        theoretical_MCSE_variance_mean_d=math.sqrt(max(0.,(mu4-((r-3)/(r-1))*v*v)/r)),
        theoretical_MCSE_event_count_variance=math.sqrt(max(0.,(count_fourth-((r-3)/(r-1))*count_variance**2)/r)))

def coverage_summary(covered,total):
    p=covered/total;z=NormalDist().inv_cdf(.975);a=z*z
    center=(p+a/(2*total))/(1+a/total)
    radius=z*math.sqrt(p*(1-p)/total+a/(4*total*total))/(1+a/total)
    return dict(proportion=p,binomial_MCSE=math.sqrt(p*(1-p)/total),
        MC_95pct_Wilson_interval=[max(0.,center-radius),min(1.,center+radius)])

def canonical(value):
    return (json.dumps(value,ensure_ascii=False,indent=2,allow_nan=False)+'\n').encode()

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def log_text(path):
    prefix=r'^\ufeff?\d{4}-\d\d-\d\dT\d\d:\d\d:\d\d(?:\.\d+)?Z ?'
    return '\n'.join(re.sub(prefix,'',line) for line in Path(path).read_text().splitlines())+'\n'

def marked(log,name):
    start='REP34_FILE_BEGIN '+name+'\n';end='REP34_FILE_END '+name
    if log.count(start)!=1 or log.count(end)!=1:raise ValueError('Nonunique/incomplete34 marked block '+name)
    return log.split(start,1)[1].split(end,1)[0]

def subset_equal(expected,actual):
    if isinstance(expected,dict):
        return isinstance(actual,dict) and all(k in actual and subset_equal(v,actual[k]) for k,v in expected.items())
    if isinstance(expected,list):
        return isinstance(actual,list) and len(expected)==len(actual) and all(subset_equal(e,a) for e,a in zip(expected,actual))
    if type(expected) is float:return close(expected,actual,rel=1e-10,absolute=1e-14)
    return type(expected) is type(actual) and expected==actual

class Audit:
    def __init__(self):self.checks=[];self.summary={}
    def check(self,name,passed,detail=None):
        self.checks.append(dict(name=name,passed=bool(passed),**({'detail':detail} if detail is not None else {})))

def expected_qualification(candidate,required):
    blockers=[dict(check=name,**candidate['checks'][name]) for name in required if candidate['checks'][name]['state']!='verified_true']
    if candidate['source_kind']!='person_level_release':
        blockers.append(dict(check='source_kind',state='verified_false',source_ids=[],reason='Declared '+candidate['source_kind']+' is not a person-level release.'))
    access=candidate['access_state']
    if access!='public_anonymous':
        reasons={'controlled_not_requested':'Controlled access has not been requested or granted.',
            'provider_prohibited_ai':'Provider prohibits the proposed AI processing.',
            'terms_or_registration_required':'New terms or registration are required and have not been granted.',
            'unknown':'Anonymous public processing access has not been verified.'}
        blockers.append(dict(check='access_state',state='not_verified' if access=='unknown' else 'verified_false',source_ids=[],reason=reasons[access]))
    return blockers

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--repo',default='.')
    p.add_argument('--log',default='/tmp/rep34_success_job.log')
    p.add_argument('--out',default='/tmp/rep34_audit')
    p.add_argument('--input-dir')
    p.add_argument('--source-execution-origin',choices=('ci','local'),default='ci')
    p.add_argument('--science-commit',required=True)
    p.add_argument('--run-id',type=int,required=True)
    p.add_argument('--job-id',type=int)
    p.add_argument('--tests-per-python',type=int,required=True)
    p.add_argument('--ci-jobs',default='/tmp/rep34_independent_ci_jobs.json')
    p.add_argument('--ci-log-dir',default='/tmp')
    p.add_argument('--ci-log-prefix',default='rep34_independent_job_')
    args=p.parse_args();repo=Path(args.repo);mortality=repo/'research/mortality';root=mortality/'replication34'
    a=Audit();log=log_text(args.log)
    names=['RESULTS_SUMMARY.json','PRECISION_RESULTS.json','SIMULATION_VALIDATION.json','READINESS_RESULTS.json','SOURCE_MANIFEST.json']
    input_dir=Path(args.input_dir) if args.input_dir else root
    paths={name:input_dir/name for name in names};data={}
    for name,path in paths.items():
        text=marked(log,name);a.check(name+' exact saved/production-stdout bytes',path.read_bytes()==text.encode());data[name]=json.loads(text)
        if input_dir!=root:a.check(name+' exact repository/localproduction bytes',(root/name).read_bytes()==path.read_bytes())
    summary=data['RESULTS_SUMMARY.json'];plans=data['PRECISION_RESULTS.json'];simulations=data['SIMULATION_VALIDATION.json'];readiness=data['READINESS_RESULTS.json'];sources=data['SOURCE_MANIFEST.json']
    lock=json.loads((root/'LOCK.json').read_text());source_lock=json.loads((root/'SOURCE_LOCK.json').read_text())
    scenario_input=json.loads((root/'SCENARIOS.json').read_text());simulation_input=json.loads((root/'SIMULATION_CASES.json').read_text())
    candidates=json.loads((root/'CANDIDATES.json').read_text());gate_contract=json.loads((root/'GATE_CONTRACT.json').read_text())
    a.check('immutable protocol ancestor of scientific commit',subprocess.run(['git','merge-base','--is-ancestor',PROTOCOL_COMMIT,args.science_commit],cwd=repo,capture_output=True).returncode==0)
    a.check('summary protocol phase schema',summary['schema']=='replication34-summary/1' and summary['plan_commit']==PROTOCOL_COMMIT and summary['phase']=='replication_readiness_and_conditional_precision_only')
    for name,digest in lock['metadata_sha256'].items():
        path=root/name if name in lock['local_metadata'] else mortality/name
        a.check('metadata pin '+name,sha(path)==digest)
    old=lock['historical_scientific_python_sha256'];actual={str(x.relative_to(mortality)):sha(x) for x in sorted(mortality.rglob('*.py'))}
    a.check('103 old scientific module hashes unchanged',len(old)==summary['prior_scientific_files_verified']==103 and all(actual.get(k)==v for k,v in old.items()))
    a.check('complete107 current scientific module hashes',len(actual)==107 and summary['scientific_python_sha256']==actual)
    a.check('four additions exclusively replication34',len(set(actual)-set(old))==4 and all(k.startswith('replication34/') for k in set(actual)-set(old)))
    a.check('old33B neutral result exact preserved',sha(mortality/'cause33b/RESULTS_SUMMARY.json')=='952a23344a6d16a453654de50bc4af437da605df9f8b2d120ced58fbbc9ec855')
    a.check('current metadata source locks exact',summary['metadata_lock_sha256']==sha(root/'LOCK.json') and summary['source_lock_sha256']==sha(root/'SOURCE_LOCK.json') and
        summary['source_manifest_sha256']==sha(root/'SOURCE_MANIFEST.json') and len(lock['metadata_sha256'])==22)
    for name,digest in source_lock['files_sha256'].items():a.check('source metadata file SHA '+name,sha(root/name)==digest)
    for name,digest in summary['outputs_sha256'].items():a.check('completed result output SHA '+name,sha(paths[name])==digest)
    a.check('no real outcome model or laboratory change claims',all(summary[k] is False for k in (
        'statistical_adequacy_for_refitted_algorithm_established','guaranteed_recruitment_or_clinical_threshold_provided',
        'new_person_level_data_or_outcomes_read','model_fits_performed','split_search_performed','old33B_result_changed','lab_report_v3_or_presentation_touched')) and
        summary['joint_integration_and_Fine_Gray']=='deferred_not_implemented_or_run' and summary['before_new_outcomes_named_release_and_precision_goals_required'] is True)
    # Source classifications are independently replayed from cited input cells.
    refs={x['id']:x for x in sources['sources']}
    a.check('18 unique named release routes and32 unique source references',len(candidates)==len({x['candidate_id'] for x in candidates})==18 and len(refs)==len(sources['sources'])==32)
    a.check('source original-byte/data claims remain false',sources['original_source_bytes_saved'] is False and sources['new_outcomes_read'] is False and
        sources['person_level_downloads_performed'] is False and all(x['original_source_sha256'] is None for x in refs.values()))
    for name,digest in sources['review_sha256'].items():a.check('source review file SHA '+name,sha(root/name)==digest)
    a.check('retained30 catalogue pin',sha(root/sources['prior_catalog']['path'])==sources['prior_catalog']['sha256'])
    a.check('current two-protein algorithm contract no falsebrand/BTP/all10 requirements',gate_contract['BTP_required'] is False and
        gate_contract['same_manufacturer_required_for_algorithm_refit'] is False and gate_contract['all_ten_NCHS_cause_codes_required'] is False and
        gate_contract['numerical_precision_approval_exists'] is False and gate_contract['no_download_or_fit_actions'] is True)
    requirements=gate_contract['requirements'];transport=gate_contract['transport_additional_requirements']
    a.check('all306 evidence cells typed cited and preserved',len(requirements)==15 and len(transport)==2 and all(
        set(c['checks'])==set(requirements+transport) and all(e['state'] in ('verified_true','verified_false','not_verified') and
            isinstance(e['reason'],str) and e['reason'].strip() and len(e['source_ids'])==len(set(e['source_ids'])) and all(k in refs for k in e['source_ids']) and
            (e['state']!='verified_true' or bool(e['source_ids'])) for e in c['checks'].values()) for c in candidates))
    action={'public_anonymous':'not_requested','controlled_not_requested':'not_requested','provider_prohibited_ai':'prohibited','terms_or_registration_required':'not_granted','unknown':'not_verified'}
    grant={'public_anonymous':'not_required_for_anonymous_public_release','controlled_not_requested':'not_granted','provider_prohibited_ai':'prohibited','terms_or_registration_required':'not_granted','unknown':'not_verified'}
    ready_by_mode={}
    for mode in ('algorithm_replication','locked_source_transport'):
        out=readiness[mode];required=requirements+(transport if mode=='locked_source_transport' else []);ready=[]
        a.check(mode+' ordered complete qualification rows',out['mode']==mode and len(out['qualifications'])==out['candidate_count']==18 and
            [x['candidate_id'] for x in out['qualifications']]==[x['candidate_id'] for x in candidates])
        for c,row in zip(candidates,out['qualifications']):
            blockers=expected_qualification(c,required);is_ready=not blockers
            if is_ready:ready.append(c['candidate_id'])
            a.check(mode+' '+c['candidate_id']+' exact independentblockers/tri-state',row['required_checks']==required and row['checks']==c['checks'] and
                row['blockers']==blockers and row['release_id']==c['release_id'] and row['source_kind']==c['source_kind'] and row['access_state']==c['access_state'] and
                row['source_ready']==row['source_ready_for_protocol']==is_ready and row['user_actions']==c['user_actions'])
            a.check(mode+' '+c['candidate_id']+' no authorization from metadata',row['access_action_state']==action[c['access_state']] and row['access_grant_state']==grant[c['access_state']] and
                row['access_actions_performed'] is False and row['person_level_download_authorized'] is False and row['outcome_comparison_authorized'] is False and
                row['precision_adequacy']=='not_established' and row['algorithm_precision_unqualified']==(mode=='algorithm_replication') and
                row['protocol_activation_status']==gate_contract['protocol_activation_status'])
        a.check(mode+' qualified count IDs preserve failclosed decision',out['source_ready_candidate_ids']==ready and out['qualified_candidate_count']==len(ready) and
            out['source_ready']==bool(ready) and out['precision_adequacy']=='not_established' and out['outcome_comparison_authorized'] is False and out['person_level_download_authorized'] is False and out['access_actions_performed'] is False)
        ready_by_mode[mode]=ready
    a.check('summary source counts and blocked activation',summary['candidates_reviewed']==18 and summary['qualified_algorithm_source_ids']==ready_by_mode['algorithm_replication'] and
        summary['qualified_transport_source_ids']==ready_by_mode['locked_source_transport'] and readiness['source_qualification_is_not_statistical_adequacy'] is True and readiness['external_predictions_or_outcome_comparisons_performed'] is False)
    # Independent scenario arithmetic, including exact raw-n binomial recurrence.
    a.check('12 frozen inputs no33B estimate defaults',len(scenario_input['scenarios'])==len(plans['cases'])==summary['illustrative_precision_scenarios']==12 and
        scenario_input['not_estimated_from_33B_effects_or_SEs'] is True and plans['input_metadata']=={k:v for k,v in scenario_input.items() if k!='scenarios'} and
        plans['estimates_are_stipulated_scenarios_not_empirical_findings'] is True and plans['adequate_sample_size_or_power_established'] is False)
    plan_rows={}
    for given,row in zip(scenario_input['scenarios'],plans['cases']):
        label=given['id'];result=row['result'];expected=expected_plan(given);plan_rows[label]=result
        a.check(label+' exact locked assumptions and ID',row['id']==label and row['assumptions']==given and result['scenario_id']==label)
        numeric={k:v for k,v in expected.items() if not k.startswith('event_')}
        a.check(label+' independent normal/variance/ceiling/yield arithmetic',subset_equal(numeric,result) and
            result['n_evaluable']==given['n_evaluable'] and result['deff_d']==given['deff_d'] and result['alpha']==given['alpha'] and result['requested_half_width']==given['half_width'])
        a.check(label+' precision honest independentofeffect andavailability',result['mean_effect_used_to_choose_n'] is False and result['algorithm_precision_status']=='adequacy_unqualified' and
            result['algorithm_training_uncertainty_included'] is False and result['full_CIF_validation_eligible'] is False and result['clinical_adequacy_or_power_established'] is False and
            result['calibration_or_AUC_precision_established'] is False and result['normal_interval_is_approximation_not_coverage_guarantee'] is True and
            result['recruitment_is_expectation_not_assurance'] is True and result['recruitment_MCAR_accounting_assumption_not_selection_correction'] is True and
            result['artifact_or_assay_qualification_performed'] is False and result['frozen_source_artifact_required_before_external_outcomes']==(given['validation_mode']=='fixed_scorer'))
        a.check(label+' censoring stop notdoublecounted endpointyield',result['IPCW_method_extension_required']==(given['censoring_fraction']>0) and
            result['bounded_uncensored_loss_formula_applies_to_full_CIF']==(given['censoring_fraction']==0) and
            result['joint_evaluable_yield_probability']==given['q_joint']*given['q_known_conditional'])
        events=result['event_count'];iid=given['event_count_model']=='independent_binomial'
        a.check(label+' independent exactbinomial quantiles actualrawN',events['model']==given['event_count_model'] and events['n_raw_trials']==given['n_evaluable'] and
            close(events['expected_raw_events'],expected['event_expected']) and close(events['variance_proxy'],expected['event_variance']) and
            events['predictive_range']==expected['event_quantiles'] and events['predictive_range_exact_distribution_assumption']==iid and
            events['predictive_range_is_not_incidence_CI'] is True and events['n_over_DEFF_used_as_binomial_trials'] is False and events['technical_floor_establishes_power_or_adequacy'] is False)
        bp=expected['event_probability_below_floor'];actual_bp=events['probability_below_technical_floor']
        a.check(label+' independent eventfloor lower-tail probability',bp is None and actual_bp is None or
            bp is not None and actual_bp is not None and math.isclose(bp,actual_bp,rel_tol=1e-9,abs_tol=1e-300))
    a.check('hypotheticalbenefit/harm/null samefixedvariance precisionN',len({plan_rows[k]['conditional_n_for_half_width'] for k in ('CENTRAL_TECHNICAL_NULL','HYPOTHETICAL_BENEFIT','HYPOTHETICAL_HARM')})==1)
    a.check('fourtimesN halvesconditionalSE',close(plan_rows['LARGER_EVALUABLE_DOMAIN']['conditional_SE'],plan_rows['CENTRAL_TECHNICAL_NULL']['conditional_SE']/2))
    # Four synthetic distributions only; never rerun producer MonteCarlo arrays.
    a.check('four frozen simulationcases andreplicate budgets',len(simulation_input['cases'])==len(simulations['cases'])==summary['synthetic_validation_cases']==4 and
        simulations['input_metadata']=={k:v for k,v in simulation_input.items() if k!='cases'} and simulations['fitted_algorithm_precision_or_external_validation_performed'] is False)
    all_moment_flags=[];sim_summaries=[]
    for given,row in zip(simulation_input['cases'],simulations['cases']):
        label=given['id'];result=row['result'];expected=expected_simulation(given);replicates=given['replicates'];n=given['n_full_clusters']*given['cluster_size']
        a.check(label+' exactfrozen atoms/clusters/seed/20kreplicates',row['id']==label and row['assumptions']==given and result['seed']==given['seed'] and result['replicates']==replicates==20000 and
            result['n_full_clusters']==given['n_full_clusters'] and result['cluster_size']==given['cluster_size'] and result['rho']==given['rho'] and result['alpha']==given['alpha'])
        exact=expected['exact_moments'];reported=result['exact_moments']
        a.check(label+' independent rationalatom moments+variance decomposition',subset_equal({k:v for k,v in exact.items() if k!='fourth_central_moment_d'},reported) and
            reported['atoms_n']==len(given['atoms']) and close(reported['variance_d'],reported['mean_conditional_Y_variance']+reported['variance_of_conditional_means']) and
            close(reported['variance_decomposition_error'],abs(reported['variance_d']-reported['mean_conditional_Y_variance']-reported['variance_of_conditional_means'])) and
            reported['predictions_fixed_before_Y'] is True and reported['fitted_model_or_clinical_validation'] is False)
        theory_keys=['n_participants','deff','exact_mean_d','exact_variance_mean_d','exact_normal_half_width','MCSE_mean_d','exact_expected_event_count','exact_event_count_variance_in_this_toy','MCSE_mean_event_count']
        a.check(label+' independent clusterICC DEFF/loss+event theory',subset_equal({k:expected[k] for k in theory_keys},result) and result['n_participants']==n)
        variance=result['MC_variance_mean_d'];fourth=result['MC_centered_fourth_moment_dbar']
        mcse_var=math.sqrt(max(0.,(fourth-((replicates-3)/(replicates-1))*variance**2)/replicates))
        a.check(label+' MCvariance/fourthmoment/SE arithmetic',variance>=0 and fourth>=0 and close(result['estimated_MCSE_mean_d'],math.sqrt(variance/replicates)) and
            close(result['approximate_fourth_moment_MCSE_variance'],mcse_var) and result['MC_variance_event_count']>=0 and -1<=result['MC_mean_d']<=1 and 0<=result['MC_mean_event_count']<=n)
        for key in ['known_toy_variance_normal_CI_coverage','estimated_cluster_mean_SE_normal_CI_coverage']:
            coverage=result[key]
            if coverage is None:
                a.check(label+' '+key+' onlyundefined forsinglecluster',given['n_full_clusters']==1);continue
            count=coverage['covered'];total=coverage['replicates'];expected_cov=coverage_summary(count,total)
            a.check(label+' '+key+' independentWilson/MCSEcounts',type(count) is int and 0<=count<=total==replicates and subset_equal(expected_cov,coverage))
        checks=row['moment_checks'];moment_specs={
            'mean_D':(result['MC_mean_d'],result['exact_mean_d'],result['MCSE_mean_d']),
            'variance_Dbar':(result['MC_variance_mean_d'],result['exact_variance_mean_d'],result['approximate_fourth_moment_MCSE_variance']),
            'mean_event_count':(result['MC_mean_event_count'],result['exact_expected_event_count'],result['MCSE_mean_event_count'])}
        for key,(observed,expected_value,se) in moment_specs.items():
            standardized=None if se==0 else abs(observed-expected_value)/se
            passed=observed==expected_value if se==0 else standardized<=6
            all_moment_flags.append(passed);check=checks[key]
            a.check(label+' '+key+' predeclaredsixMCSEnumericcheck',close(check['observed'],observed) and close(check['expected'],expected_value) and close(check['MCSE'],se) and
                close(check['absolute_discrepancy_in_MCSE'],standardized) and check['predeclared_six_MCSE_numerical_check_passed']==passed)
        a.check(label+' noalgorithmclinicalclaim orpositivevariance tolerance',row['coverage_95pct_is_reported_not_required_or_tuned'] is True and
            result['positive_variance_CI_coverage_uses_no_added_tolerance'] is True and result['predictions_fixed_before_Y'] is True and
            result['refitted_algorithm_precision_validated'] is False and result['external_cohort_or_clinical_validation_performed'] is False and result['participant_simulation_arrays_exported'] is False and
            result['event_and_loss_DEFF_equal_only_under_this_specified_shared_pair_law'] is True)
        if expected['exact_variance_mean_d']==0:
            a.check(label+' deterministiczero-loss case not95coveragevalidation',result['MC_mean_d']==result['MC_variance_mean_d']==0 and result['known_toy_variance_normal_CI_coverage']['covered']==replicates)
        sim_summaries.append(dict(id=label,exact_mean=expected['exact_mean_d'],exact_variance_mean=expected['exact_variance_mean_d'],
            MC_mean=result['MC_mean_d'],MC_variance_mean=variance,known_coverage=result['known_toy_variance_normal_CI_coverage'],estimated_coverage=result['estimated_cluster_mean_SE_normal_CI_coverage']))
    a.check('summary all12numeric moment flags reported withoutcoveragecriterion',summary['moment_checks_passed']==simulations['all_predeclared_moment_checks_passed']==all(all_moment_flags) and
        summary['MonteCarlo_replicates']=={x['id']:x['replicates'] for x in simulation_input['cases']})
    # Final CI evidence is separate from production origin; push-only jobs do
    # not turn a safe local scenario run into a scenario-CI success.
    if Path(args.ci_jobs).exists():
        jobs=json.loads(Path(args.ci_jobs).read_text())
        unit_jobs=[j for j in jobs if j['name'].startswith('unit')]
        scenario_jobs=[j for j in jobs if not j['name'].startswith('unit')]
        a.check('completed scientificCI unitjobs success',len(unit_jobs)==2 and all(j['run_id']==args.run_id and j['status']=='completed' and j['conclusion']=='success' for j in unit_jobs))
        if args.source_execution_origin=='local':
            a.check('localproduction no false scenarioCI success',len(scenario_jobs)==1 and all(j['status']=='completed' and j['conclusion']=='skipped' for j in scenario_jobs))
        else:
            a.check('scenarioCI production success',len(scenario_jobs)==1 and all(j['run_id']==args.run_id and j['id']==args.job_id and j['status']=='completed' and j['conclusion']=='success' for j in scenario_jobs))
        for job in jobs:
            if job['conclusion']=='skipped':continue
            ownlog=log_text(Path(args.ci_log_dir)/(args.ci_log_prefix+str(job['id'])+'.log'))
            a.check(job['name']+' scientifichead incheckout',args.science_commit in ownlog)
            if job['name'].startswith('unit'):
                a.check(job['name']+' syntheticunitcount',bool(re.search(r'Ran '+str(args.tests_per_python)+r' tests in [\d.]+s\n\s*\nOK\n',ownlog)))
            else:a.check('five exports match independentcompletedscenarioslog',all(marked(ownlog,name)==marked(log,name) for name in names))
    else:a.check('independent completedCI evidence supplied',False)
    def forbidden(value):
        if isinstance(value,dict):return [k for k in value if k.lower() in {'seqn','participant_ids','individual_labels','individual_predictions','participant_predictions','oof_predictions','simulated_participant_arrays'}]+[k for v in value.values() for k in forbidden(v)]
        if isinstance(value,list):return [k for v in value for k in forbidden(v)]
        return []
    a.check('onlymetadata/scenarios/aggregate simulations exported',not any(forbidden(x) for x in data.values()))
    merge=re.search(r'git log -1 --format=%H\n([0-9a-f]{40})',log)
    a.summary=dict(protocol_commit=PROTOCOL_COMMIT,science_commit=args.science_commit,run_id=args.run_id,job_id=args.job_id,
        source_execution_origin=args.source_execution_origin,
        actual_PR_merge_checkout=merge.group(1) if merge and args.source_execution_origin=='ci' else None,
        scenario_CI_executed=args.source_execution_origin=='ci',
        local_exit0_reported_by_parent=args.source_execution_origin=='local',tests_per_python=args.tests_per_python,
        candidates18=len(candidates),sources32=len(refs),source_ready_ids=ready_by_mode,
        evidence_state_counts=dict(Counter(e['state'] for c in candidates for e in c['checks'].values())),
        scenario_count=len(plans['cases']),simulation_count=len(simulations['cases']),all_predeclared_moment_checks_passed=all(all_moment_flags),
        central_technical_plan=plan_rows['CENTRAL_TECHNICAL_NULL'],scenario_plans=plan_rows,simulation_summaries=sim_summaries,
        source_files_sha256={k:sha(root/k) for k in source_lock['files_sha256']},output_sha256={k:sha(v) for k,v in paths.items()},
        scientific_files=len(actual),old_scientific_files=len(old),production_stdout_sha256=sha(args.log))
    limits=[
        'No new cohort records, outcomes, empirical model fits, old33B retuning or external predictions were read or performed. Qualification is metadata-only; all current routes remain source-blocked or unverified.',
        'Conditional paired-loss precision is not adequacy for a refitted algorithm. All scenario effects/SD/SE/incidences/halfwidths are predeclared technical assumptions, not estimates from33B, clinical thresholds, guaranteedN or posthocpower.',
        'Exact binomial predictive counts use actualraw evaluable n only under independentBernoulli assumption. Weighted national incidence and lossDEFF do not define raw sample-event distributions or AUC/calibration precision.',
        'Availability is jointpanel/core and endpointknown conditionalyield; recruitment is expectation-only. Anybefore96censoring requires a separateIPCW method extension and does not become valid by recruiting more or dropping censored people.',
        'Synthetic finitepredictionatoms/cluster-mixture validate only mathematical implementation under thatspecifiedDGP. Reported normalCI coverage, includingrare/degenerate cases, is not forcedto95 and is not clinical/algorithm validation.',
        'The auditor independently recomputes scenarios, relative-PMFbinomial quantiles, exact rationalatom moments, clusterDEFF,theoreticalmoments,Wilsoncoverage/MCSE and allqualificationdecisions. It doesnot rerun the80,000 MonteCarlo datasets or reconstruct their original random draws.',
        'Sourceclassifications useprovidedprimary-source reviews and citedmetadata; externaloriginalsource bytes were not independentlydownloaded/rehashed. Savedreview/ledger hashes and exactproduction plaintext exports are checked; noZIP/raw retrieval.',
        'Source-ready wouldmean readiness for a laternamedrelease protocol, not authorizeddownload/comparison, sufficientpower, clinicalreadiness or a frozen-source artifact. Oldmodels, neutral33B, labreportv3/presentation anddeferredFG/jointscope are preserved.'
    ]
    if args.source_execution_origin=='local':
        limits.append('Actual5JSON production was executed locally byroot at the scientific commit and audited from its savedstdout/files; exit0 is reported byroot. PushCI independently confirms unit41 tests perversion only; its PR-only scenarios job wasSKIPPED. No scenarioCI pass,PRmerge checkout or independent rerun is claimed. DraftPR creation was rejected byautomaticapprovalreview; no retry/workaround was used.')
    report=dict(schema='independent-replication34-audit/1',passed=sum(x['passed'] for x in a.checks),total=len(a.checks),all_passed=all(x['passed'] for x in a.checks),summary=a.summary,checks=a.checks,limits=limits)
    out=Path(args.out);out.with_suffix('.json').write_bytes(canonical(report))
    lines=['# Независимый аудит34','',f"**{report['passed']}/{report['total']} PASS; общий исход: {'PASS' if report['all_passed'] else 'FAIL'}.**",'',
        f"Protocol `{PROTOCOL_COMMIT}`; scientific commit `{args.science_commit}`; production origin **{args.source_execution_origin}**; unitCIrun{args.run_id}.",
        (f"ScenariosCIjob{args.job_id}; actual PRmerge `{a.summary['actual_PR_merge_checkout']}`." if args.source_execution_origin=='ci' else
         'ПятьJSON созданы локальным запускомroot; exactstdout/localfiles/repositorybytes проверены. ScenariosCIjob SKIPPED, PRmerge не заявлен; повторного запускаMC аудитором не было.'),
        f"CI: {args.tests_per_python} unit tests on each Python3.11/3.13.18routes/32references; ready algorithm={ready_by_mode['algorithm_replication']},transport={ready_by_mode['locked_source_transport']}.",'',
        'Все12планов пересчитаны независимой арифметикой. Четыре toyDGP проверены exactrational moments,clusterDEFF иMonteCarlo/Wilson identities; simulations не повторялись.','',
        '| Scenario | Conditional SE | Conditional halfwidth | Conditional n target | Expected-yield recruitment |','|---|---:|---:|---:|---:|']
    for label,result in plan_rows.items():lines.append(f"| {label} | {result['conditional_SE']} | {result['conditional_normal_half_width']} | {result['conditional_n_for_half_width']} | {result['candidate_recruitment_expectation']} |")
    lines+=['','ВсеN — conditionalfixed-loss technicalproxies, **не adequateN дляalgorithm replication**, неpower/гарантияevents инеclinicalthresholds.','',
        '| Synthetic case | Exact mean D | Exact variance mean D | MC mean D | MC variance | Known-variance coverage | Estimated-cluster-SE coverage |',
        '|---|---:|---:|---:|---:|---:|---:|']
    for row in sim_summaries:
        ec=row['estimated_coverage'];lines.append(f"| {row['id']} | {row['exact_mean']} | {row['exact_variance_mean']} | {row['MC_mean']} | {row['MC_variance_mean']} | {row['known_coverage']['proportion']} | {ec['proportion'] if ec else None} |")
    lines+=['','Coverage тольковымышленногозакона; rare/degeneratecase не обязан иметь95%. Sourceblockers остаютсяотдельно отsuccessfulsyntheticmath.','',
        'Ограничения:','']+['- '+x for x in limits]
    failed=[x['name'] for x in a.checks if not x['passed']]
    if failed:lines+=['','Failed checks:','']+['- '+x for x in failed]
    out.with_suffix('.md').write_text('\n'.join(lines)+'\n')
    print(json.dumps(dict(passed=report['passed'],total=report['total'],all_passed=report['all_passed'],failed=failed),ensure_ascii=False))
    return 0 if report['all_passed'] else 1

if __name__=='__main__':raise SystemExit(main())
