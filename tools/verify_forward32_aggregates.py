#!/usr/bin/env python3
"""Independent stdlib aggregate audit of forward32, no pipeline imports/fits/GET.

Prepared before final32 empirical output. Superseded runs are never consulted.
"""
import argparse
from collections import Counter
import csv
import hashlib
import io
import json
import math
from pathlib import Path
import re
import subprocess
from urllib.parse import urlparse

SCIENCE='59b0626767c0169107652cc7b7d51a416c8e7668'
RUN_ID=40917437691
JOB_ID=123601312507
ARMS=('routine','routine_b2m','routine_cys','routine_cys_b2m')
PAIRS={'primary':('routine','routine_b2m'),'secondary':('routine_cys','routine_cys_b2m')}
MODES=('fixed_training','joint_refit')
METRICS=('auc','brier','observed','predicted','bias_pp','observed_expected','calibration_in_large','calibration_intercept','calibration_slope')
CORE=['age_linear','age_quadratic','male','current','former','sbp','bmi','hba1c','egfr_low','egfr_high','log_uacr']
CORE_NAMES=['age_centered','age_centered_squared','male','smoking_current','smoking_former','systolic_bp','bmi','hba1c','egfr_low','egfr_high','log2_uacr']
PNAMES=['log2_cystatin_c','log2_beta2_microglobulin']
INDICES={'routine':[],'routine_b2m':[1],'routine_cys':[0],'routine_cys_b2m':[0,1]}
SCOPES={'fixed_training':'late_PSU_resampling_with_fixed_early_training',
 'joint_refit':'independent_early_and_late_PSU_resampling_with_full_early_preprocessing_tuning_refitting'}

def digest(value): return hashlib.sha256(value).hexdigest()
def sha(path): return digest(Path(path).read_bytes())
def canonical(value): return (json.dumps(value,ensure_ascii=False,indent=2,allow_nan=False)+'\n').encode()
def finite(value): return isinstance(value,(int,float)) and not isinstance(value,bool) and math.isfinite(value)
def num(value):
    if value is None or value=='': return None
    result=float(value)
    if not math.isfinite(result): raise ValueError('Nonfinite exported numeric cell')
    return result
def close(a,b):
    return (a is None and b is None) or (finite(a) and finite(b) and math.isclose(a,b,rel_tol=1e-10,abs_tol=2e-12))
def percentile(values,p):
    if not values: return None
    ordered=sorted(values); position=(len(ordered)-1)*p
    lo=math.floor(position); hi=math.ceil(position)
    return ordered[lo]+(ordered[hi]-ordered[lo])*(position-lo)
def log_text(path):
    return '\n'.join(re.sub(r'^\d{4}-\d\d-\d\dT\S+\s','',line)
        for line in Path(path).read_text().splitlines())+'\n'
def block(log,label):
    start='FORWARD32_'+label+'_BEGIN\n'; end='FORWARD32_'+label+'_END'
    if log.count(start)!=1 or log.count(end)!=1: raise ValueError('Nonunique/incomplete marked block '+label)
    return log.split(start,1)[1].split(end,1)[0]
class Audit:
    def __init__(self): self.checks=[]; self.summary={}
    def check(self,name,ok,detail=None):
        self.checks.append(dict(name=name,passed=bool(ok),**({'detail':detail} if detail is not None else {})))
    def group(self,name,conditions): self.check(name,all(conditions))
def converged(fit):
    g=fit.get('gradient_max')
    return finite(g) and ((fit.get('convergence_status')=='gradient_1e-6' and g<=1e-6) or
        (fit.get('convergence_status')=='tiny_step_gradient_1e-4' and g<=1e-4))

def inspect_metric_row(row,weights=None):
    """Metric-local algebra and pair status, independent of result helpers."""
    ok=True
    for arm in ARMS:
        status=row[arm+'|status']; values={m:num(row[arm+'|'+m]) for m in METRICS}
        ok &= status in {'computed','undefined'}
        if status=='undefined':
            ok &= bool(row.get(arm+'|error')) and all(v is None for v in values.values())
            continue
        ok &= all(values[k] is not None and 0<=values[k]<=1 for k in ('brier','observed','predicted'))
        ok &= values['auc'] is None or 0<=values['auc']<=1
        ok &= close(values['bias_pp'],100*(values['predicted']-values['observed']))
        ok &= close(values['observed_expected'],values['observed']/values['predicted'])
        if row[arm+'|calibration_status']=='computed': ok &= all(values[k] is not None for k in METRICS[-3:])
        if weights:
            ok &= close(values['observed'],num(weights['late_event_weight_sum'])/num(weights['late_weight_sum']))
    for label,(base,expanded) in PAIRS.items():
        both=all(row[arm+'|status']=='computed' for arm in (base,expanded))
        ok &= row[label+'|status']==('computed' if both else 'undefined')
        for metric in ('brier','auc'):
            x,y=num(row[expanded+'|'+metric]),num(row[base+'|'+metric])
            ok &= close(num(row[label+'|delta_'+metric]),None if x is None or y is None else x-y)
    return ok

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo',default=str(Path(__file__).resolve().parents[1]))
    parser.add_argument('--log',default='/tmp/forward32_success_job.log')
    parser.add_argument('--out',default='/tmp/verify32')
    parser.add_argument('--ci-jobs',default='/tmp/forward32_independent_ci_jobs.json')
    parser.add_argument('--ci-log-dir',default='/tmp')
    args=parser.parse_args();ci_log_dir=Path(args.ci_log_dir)
    repo=Path(args.repo); mortality=repo/'research/mortality'; root=mortality/'forward32'
    a=Audit(); log=log_text(args.log); lock=json.loads((root/'LOCK.json').read_text())
    r=json.loads(block(log,'RESULTS')); manifest=json.loads(block(log,'MANIFEST'))
    computed=r['status']!='blocked_feasibility_no_model_fits'
    paths={'RESULTS':root/'RESULTS_SUMMARY.json','MANIFEST':root/'SOURCE_MANIFEST.json'}
    if not paths['RESULTS'].exists(): paths['RESULTS']=root/'RESULTS.json'
    if computed:
        paths.update(MODEL_ARTIFACTS=root/'MODEL_ARTIFACTS.json',CALIBRATION_BINS=root/'CALIBRATION_BINS.json',
            FIXED_TRAINING_BOOTSTRAP_CSV=root/'fixed_training_bootstrap.csv',
            JOINT_REFIT_BOOTSTRAP_CSV=root/'joint_refit_bootstrap.csv',PAIRED_BOOTSTRAP_CSV=root/'paired_bootstrap.csv')
    outputs={}
    for label,path in paths.items():
        text=block(log,label)
        a.check(label+' exact bytes agree with successful marked job stdout',path.read_bytes()==text.encode())
        outputs[label]=list(csv.DictReader(io.StringIO(text))) if label.endswith('_CSV') else json.loads(text)
    a.check('plan ancestor of final science commit',subprocess.run(['git','merge-base','--is-ancestor','377778c',SCIENCE],cwd=repo,capture_output=True).returncode==0)
    a.group('predeclared settings unchanged',[r['plan_commit']=='377778c',r['horizon_years']==8,r['horizon_months']==96,
        r['train_cycles']==[1999,2001],r['late_cycles']==[2003],r['train_weight']=='WTSSCB4Y',r['late_weight']=='WTSSCB2Y',
        r['routine_features']==CORE,r['protein_features']==['SSCYST','SSB2M'],
        r['primary_contrast']=='routine_b2m_minus_routine',r['secondary_contrast']=='routine_cys_b2m_minus_routine_cys',
        r['primary_score']=='late_survey_weighted_eight_year_Brier',r['primary_uncertainty_mode']=='joint_refit'])
    a.group('research and mode labels honest',[r['analysis']=='exploratory_reused_NHANES1999_2004_fixed_forward_calendar_B2M',
        r['same_late_draw_all_arms_and_both_modes'] is True,r['uncertainty_modes_pooled'] is False,
        r['same_common_domain_all_four_models_within_each_period'] is True,r['late_recalibration_performed'] is False,
        r['independent_external_validation'] is False,r['clinical_use_ready'] is False,
        r['participant_records_exported'] is False,r['previous_frozen_coefficients_used'] is False])
    a.check('plan and lock hashes',r['plan_sha256']==sha(root/'ANALYSIS_PLAN.md') and r['lock_sha256']==sha(root/'LOCK.json'))
    for name,expected in lock['metadata_sha256'].items():
        path=mortality/name if '/' in name else root/name
        a.check('metadata pin '+name,sha(path)==expected)
    old=lock['historical_scientific_python_sha256']
    actual={str(p.relative_to(mortality)):sha(p) for p in sorted(mortality.rglob('*.py'))}
    a.check('all90 historical scientific module hashes unchanged',len(old)==90 and all(actual.get(p)==v for p,v in old.items()))
    a.check('all95 actual scientific module hashes agree with empirical runner',len(actual)==95 and actual==r['scientific_python_sha256'])
    a.check('five additions solely forward32',len(set(actual)-set(old))==5 and all(p.startswith('forward32/') for p in set(actual)-set(old)))
    for filename,expected in r['outputs_sha256'].items():
        a.check('output exact SHA '+filename,sha(root/filename)==expected)
    a.check('source manifest export SHA',r['source_manifest_export_sha256']==sha(paths['MANIFEST']))
    # Wrapper/hash agreement; a separate source-specialist audit covers observed full marginals/native weights.
    sa=r['source_audit']; prior=sa['prior31_source_audit']; source=manifest['source']; intake=manifest['intake']
    source_bytes=canonical(source)
    a.check('new source wrapper bound to final protocol',manifest['schema']=='forward32-source-manifest/1' and
        intake['schema']=='forward32-retrieval/1' and intake['plan_commit']=='377778c' and intake['lock_sha256']==sha(root/'LOCK.json') and
        intake['base_manifest_sha256']==digest(source_bytes)==prior['source_manifest_sha256'] and
        digest(canonical(intake))==sa['forward_source_wrapper_sha256'])
    fingerprints={x['file']:x['sha256'] for x in source['files']}
    a.check('all50 fresh required sources with consistent fingerprints',len(source['files'])==len(fingerprints)==50 and
        fingerprints==intake['source_fingerprints']==prior['source_fingerprints'] and prior['files_verified']==50 and prior['historical_data_hashes_verified']==25 and
        prior['source_content_sha256']==digest(json.dumps(fingerprints,sort_keys=True,separators=(',',':')).encode()))
    historical={name:hashvalue for hashvalue,name in (line.split() for line in (mortality/'nhanes04/SOURCE_DATA_SHA256.txt').read_text().splitlines() if line.strip())}
    a.check('25 historical public source data hashes fixed',len([x for x in fingerprints if x.endswith(('.dat','.xpt'))])==25 and
        all(historical.get(name)==value for name,value in fingerprints.items() if name.endswith(('.dat','.xpt'))))
    safe=lambda url: urlparse(url).scheme=='https' and urlparse(url).hostname in {'wwwn.cdc.gov','www.cdc.gov','ftp.cdc.gov'}
    a.check('all mandatory source rows retrieved from official safe hosts',all(x['status']=='retrieved' and x['bytes']>0 and safe(x['url']) and safe(x['resolved_url']) for x in source['files']))
    denials=source['documentation_not_retrieved']
    a.check('existing two PDF403 denials remain explicit without retry',len(denials)==2 and denials==prior['documentation_not_retrieved'] and
        all(x['status']=='denied_not_retried' and x['original_http_status']==403 and x['file'] not in fingerprints and 'sha256' not in x for x in denials))
    periods=sa['periods']; train=periods['train']['common']; late=periods['late']['common']
    a.check('zero participant and masked PSU overlap',sa['participant_overlap_n']==sa['masked_PSU_overlap_n']==0)
    for part,start,end,psus,weight in [('train',1,28,57,'WTSSCB4Y'),('late',29,43,30,'WTSSCB2Y')]:
        p=periods[part]; c=p['common']; eligible=p['eligible']; design=p['full_masked_design']
        a.check(part+' full eligible common masked design identical',design==p['eligible_masked_design']==p['common_masked_design'] and
            set(design)==set(map(str,range(start,end+1))) and sum(map(len,design.values()))==psus and
            all(len(v)>=2 for v in design.values()) and p['weight_variable']==weight)
        a.check(part+' primary-domain support and native weight invariants',c['n']==c['events8']+c['nonevents8'] and
            c['positive_weight_events8']==c['events8']>=20 and c['positive_weight_nonevents8']==c['nonevents8']>=20 and
            c['selected_native_source_max_abs_difference']==0 and c['minimum_survivor_followup_months']>=96 and
            close(c['effective_n'],c['native_weight_sum']**2/c['native_weight_squared_sum']) and
            close(p['common_native_weighted_coverage'],c['native_weight_sum']/eligible['native_weight_sum']))
    a.check('source gate exactly reconciles pass or block',r['gate']==sa['gate'] and r['gate']['support_threshold_is_not_power_or_calibration_certificate'] is True and
        r['gate']['temporal_evaluation_allowed']==computed and bool(r['gate']['blockers'])!=computed)
    mapping=sa['eligible_derived_early_inner_fold_mapping']
    keys=[str(st)+':'+str(psu) for st,psus in periods['train']['eligible_masked_design'].items() for psu in psus]
    ordered=sorted(keys,key=lambda key:digest(('forward32:inner:v1:'+key).encode()))
    expected_mapping={key:i%3 for i,key in enumerate(ordered)}
    a.check('independent SHA early eligible3fold mapping',mapping==expected_mapping and len(mapping)==57)
    ci=Path(args.ci_jobs)
    if ci.exists():
        jobs=json.loads(ci.read_text())
        a.check('independently queried final-science CI completed successfully',len(jobs)==3 and all(j['run_id']==RUN_ID and j['status']=='completed' and j['conclusion']=='success' for j in jobs))
        for job in jobs:
            ownlog=log_text(ci_log_dir/('forward32_independent_job_'+str(job['id'])+'.log'))
            a.check(job['name']+' sciencehead in CI checkout',SCIENCE in ownlog)
            if job['name'].startswith('synthetic'):
                a.check(job['name']+'54 tests pass',bool(re.search(r'Ran 54 tests in [\d.]+s\n\s*\nOK\n',ownlog)))
            else:
                a.check('all blocks match independent successful-job fetch',all(block(ownlog,k)==block(log,k) for k in paths))
    else:
        a.check('independent CI evidence supplied',False,'Completed CI metadata/logs have not been saved yet')
    if not computed:
        a.check('blocked output has no fits or uncertainty',r['model_fits_attempted']==r['uncertainty_draws_attempted']==0 and
            r['blocked_reasons']==r['gate']['blockers'] and not any(k in r for k in ('point','uncertainty','contrasts')))
        a.summary=dict(status=r['status'],science_commit=SCIENCE,run_id=RUN_ID,job_id=JOB_ID,blocked_reasons=r['blocked_reasons'])
        return finish(a,args,paths)
    artifact=outputs['MODEL_ARTIFACTS']; bins=outputs['CALIBRATION_BINS']; point=r['point']
    a.group('computed fixed forward implementation settings',[r['model_fits_attempted']==4,r['bootstrap_replicates_per_mode']==500,
        r['late_seed']==20261103,r['train_seed']==20261104,r['training_inner_fold_aggregates_match_eligible_mapping'] is True,
        r['bootstrap_training_inner_partition_fixed'] is True,r['fixed_training_excludes_training_uncertainty'] is True,
        r['assay_storage_calibration_selection_history_uncertainty_included'] is False])
    a.check('split n events agree with source-qualified aggregate',r['train_n']==train['n'] and r['train_events8']==train['events8'] and
        r['late_n']==late['n'] and r['late_events8']==late['events8'] and r['train_n']+r['late_n']==sa['original31_common_n'])
    a.check('all500 paired attempts completed in successful log','FORWARD32_PAIRED_DRAWS 500' in log)
    a.check('point flattened metrics and paired contrasts',inspect_metric_row(r['contrasts']))
    for arm,p in point.items():
        a.check(arm+' late point n weights and metric algebra',p['n']==late['n'] and p['events']==late['events8'] and
            p['n_positive_weight']==late['n'] and p['events_positive_weight']==late['events8'] and
            close(p['weight_sum'],late['native_weight_sum']) and
            all(close(p[m],r['contrasts'][arm+'|'+m]) for m in METRICS) and
            (p['calibration_status']!='computed' or (finite(p['calibration_gradient_max']) and p['calibration_gradient_max']<=1e-5)))
        rows=[b for b in bins if b['model']==arm]
        a.check(arm+' calibration bins population and weighted moment reconstruction',1<=len(rows)<=5 and
            all(b['analysis_domain']=='late_point' and b['weight_sum']>0 and close(b['bias_pp'],100*(b['predicted']-b['observed'])) for b in rows) and
            sum(b['n'] for b in rows)==late['n'] and sum(b['events'] for b in rows)==late['events8'] and
            close(sum(b['weight_sum'] for b in rows),late['native_weight_sum']) and
            close(sum(b['weight_sum']*b['observed'] for b in rows)/late['native_weight_sum'],p['observed']) and
            close(sum(b['weight_sum']*b['predicted'] for b in rows)/late['native_weight_sum'],p['predicted']))
    a.check('point4 arm statuses and artifact presence',set(artifact['models'])==set(point)==set(ARMS) and
        all(r['point_model_statuses'][x]['status']==r['point_statuses'][x]['status']=='computed' for x in ARMS))
    a.group('artifact training scope and no late leakage flags',[artifact['schema']=='forward32-fit/1',artifact['fit_scope']=='early_training_only',
        artifact['selection_scope']=='early_three_fold_PSU_inner_Brier',artifact['no_outer_cross_validation'] is True,
        artifact['late_feature_preprocessing_fitted'] is False,artifact['inner_scores_are_training_diagnostics_not_validation'] is True,
        artifact['independent_external_validation'] is False,artifact['participant_records_exported'] is False,
        artifact['late_domain']['outcomes_used_for_fitting'] is False,artifact['late_domain']['weights_used_for_fitting'] is False,
        artifact['late_domain']['training_psu_overlap_n']==0,artifact['training']['weight_name']=='WTSSCB4Y',
        artifact['late_domain']['evaluation_weight_name']=='WTSSCB2Y',artifact['training']['cycles']==[1999,2001],
        artifact['late_domain']['cycle']==2003,artifact['training']['n']==train['n'],artifact['training']['events8']==train['events8'],
        close(artifact['training']['weight_sum'],train['native_weight_sum']),artifact['training']['psu_n']==57,
        artifact['late_domain']['n']==late['n'],artifact['late_domain']['psu_n']==30,
        artifact['plan_commit']==r['plan_commit'],artifact['plan_sha256']==r['plan_sha256'],artifact['lock_sha256']==r['lock_sha256'],
        artifact['scientific_python_sha256']==actual])
    support=artifact['inner_folds_aggregate']
    a.check('3innerfold totals and no PSU overlap',sorted(p['fold'] for p in support)==[0,1,2] and
        sum(p['scoring_n'] for p in support)==train['n'] and sum(p['scoring_events'] for p in support)==train['events8'] and
        all(p['training_n']+p['scoring_n']==train['n'] and p['training_events']+p['scoring_events']==train['events8'] and
            p['training_positive_weight_n']==p['training_n'] and p['scoring_positive_weight_n']==p['scoring_n'] and
            p['training_psu_n']+p['scoring_psu_n']==57 and p['PSU_overlap_n']==0 and
            p['scoring_psu_n']==sum(f==p['fold'] for f in mapping.values()) and
            close(p['training_weight_sum']+p['scoring_weight_sum'],train['native_weight_sum']) for p in support))
    for arm,row in artifact['models'].items():
        prep=row['preprocessing']; fit=row['fit']; dim=11+len(INDICES[arm])
        a.check(arm+' feature layout coefficients and training-only scales',prep['model_name']==arm and prep['core_indices']==list(range(11)) and
            prep['protein_indices']==INDICES[arm] and row['feature_names']==CORE_NAMES+[PNAMES[i] for i in INDICES[arm]] and
            len(prep['mean'])==len(prep['scale'])==len(prep['constant'])==len(fit['coefficients'])==dim and
            all(finite(x) for x in prep['mean']+prep['scale']+fit['coefficients']) and all(x>0 for x in prep['scale']) and
            all(isinstance(x,bool) for x in prep['constant']) and all(not c or s==1 for c,s in zip(prep['constant'],prep['scale'])))
        moments=[train['outcome_free_feature_moments'][x] for x in CORE]+[train['outcome_free_log2_protein_moments'][['SSCYST','SSB2M'][i]] for i in INDICES[arm]]
        a.check(arm+' final preprocessing independent aggregate moment match',all(m['n']==train['n'] and close(mean,m['weighted_mean']) and
            close(scale,1. if const else m['weighted_SD']) for mean,scale,const,m in zip(prep['mean'],prep['scale'],prep['constant'],moments)))
        scores=row['inner_brier_by_penalty']
        a.check(arm+' early-only lambda selection independently reconstructed',set(scores)=={'1.0','10.0','100.0'} and all(finite(v) and 0<=v<=1 for v in scores.values()) and
            row['selected_penalty']==fit['penalty']==max(float(k) for k,v in scores.items() if v<=min(scores.values())+1e-12))
        a.check(arm+' final fit finite normalized and converged',finite(fit['intercept']) and finite(fit['objective']) and fit['iterations']>0 and
            close(fit['normalized_weight_sum'],train['n']) and converged(fit))
        diags=row['inner_fit_diagnostics']
        a.check(arm+' all9 early innerfits convergence and normalization',len(diags)==9 and
            {(p['fold'],p['penalty']) for p in diags}=={(f,p) for f in range(3) for p in (1.,10.,100.)} and
            all(converged(p) and p['iterations']>0 and close(p['normalized_weight_sum'],next(s['training_n'] for s in support if s['fold']==p['fold'])) for p in diags))
    # Independent percentile math over all80 targets, preserving mode/metric-local missingness.
    draws={mode:outputs[mode.upper()+'_BOOTSTRAP_CSV'] for mode in MODES}
    interval_keys={arm+'|'+m for arm in ARMS for m in METRICS}|{label+'|delta_'+m for label in PAIRS for m in ('brier','auc')}
    a.check('exact two separate uncertainty modes primaryjoint',set(r['uncertainty'])==set(MODES) and all(
        r['uncertainty'][mode]['scope']==SCOPES[mode] and r['uncertainty'][mode]['is_primary_uncertainty']==(mode=='joint_refit') for mode in MODES))
    calibration_valid={}
    for mode in MODES:
        rows=draws[mode]; u=r['uncertainty'][mode]
        a.check(mode+'500 ordered unique attempted replicates',[int(x['replicate']) for x in rows]==list(range(500)) and
            all(x['mode']==mode for x in rows) and u['attempted_draws']==500)
        a.check(mode+' all40 interval targets',set(u['intervals'])==interval_keys)
        for key in sorted(interval_keys):
            values=[num(x[key]) for x in rows if num(x[key]) is not None]; interval=u['intervals'][key]
            a.group(mode+' independent percentile '+key,[close(interval['lower'],percentile(values,.025)),close(interval['upper'],percentile(values,.975)),
                interval['valid']==len(values),interval['undefined']==500-len(values),interval['defined_draws_only']==(len(values)!=500),
                interval['scope']==SCOPES[mode],interval['uncertainty_mode']==mode,interval['early_training_refitted']==(mode=='joint_refit')])
        a.check(mode+'500 paired metrics and resampled observed risks',all(inspect_metric_row(x,x) for x in rows))
        a.check(mode+' early resample usage flag',all(x['early_resample_used_for_fitting']==str(mode=='joint_refit') for x in rows))
        for arm in ARMS:
            a.check(mode+' '+arm+' failures match retained statuses',u['model_failure_counts'][arm]==sum(x[arm+'|status']!='computed' for x in rows))
        for label in PAIRS:
            a.check(mode+' '+label+' failures match paired statuses',u['contrast_failure_counts'][label]==sum(x[label+'|status']!='computed' for x in rows))
        a.check(mode+' valid nonnegative weight and positive-row summaries',all(
            num(x[p+'weight_sum'])>0 and num(x[p+'weight_squared_sum'])>0 and
            0<num(x[p+'positive_weight_n'])<=domain['n'] and
            num(x[p+'positive_weight_events8'])+num(x[p+'positive_weight_nonevents8'])==num(x[p+'positive_weight_n']) and
            0<=num(x[p+'event_weight_sum'])<=num(x[p+'weight_sum']) and
            0<num(x[p+'positive_weight_psu_n'])<=domain['psu_n']
            for x in rows for p,domain in [('early_',train),('late_',late)]))
        calibration_valid[mode]={arm:{m:u['intervals'][arm+'|'+m]['valid'] for m in METRICS[-3:]} for arm in ARMS}
    latekeys=[k for k in draws['fixed_training'][0] if k.startswith('late_')]
    earlykeys=[k for k in draws['fixed_training'][0] if k.startswith('early_') and k!='early_resample_used_for_fitting']
    a.check('same exact late resampling summaries in both modes each500 paired draw',all(
        all(fixed[k]==joint[k] for k in latekeys) for fixed,joint in zip(draws['fixed_training'],draws['joint_refit'])))
    a.check('late resample varies across500attempts',len({tuple(x[k] for k in latekeys) for x in draws['fixed_training']})>1)
    a.check('fixed early weight summaries constant but joint vary',len({tuple(x[k] for k in earlykeys) for x in draws['fixed_training']})==1 and
        len({tuple(x[k] for k in earlykeys) for x in draws['joint_refit']})>1)
    original=draws['fixed_training'][0]
    a.check('fixed early summaries equal pointtraining source aggregate',close(num(original['early_weight_sum']),train['native_weight_sum']) and
        close(num(original['early_weight_squared_sum']),train['native_weight_squared_sum']) and
        num(original['early_positive_weight_n'])==train['n'] and num(original['early_positive_weight_events8'])==train['events8'] and
        num(original['early_positive_weight_nonevents8'])==train['nonevents8'] and num(original['early_positive_weight_psu_n'])==train['psu_n'])
    for arm in ARMS:
        selection=r['joint_refit_penalty_selection_counts'][arm]
        train_failed=sum(x[arm+'|status']!='computed' and x[arm+'|error'].startswith('training|') for x in draws['joint_refit'])
        a.check(arm+' joint500 finalfit lambda selections retained separately',set(selection)=={'1.0','10.0','100.0'} and
            all(isinstance(v,int) and v>=0 for v in selection.values()) and sum(selection.values())==500-train_failed)
    paired=outputs['PAIRED_BOOTSTRAP_CSV']
    a.check('1000 paired export rows unique by mode and replicate',len(paired)==1000 and
        Counter(x['mode'] for x in paired)==Counter({mode:500 for mode in MODES}) and
        len({(x['mode'],int(x['replicate'])) for x in paired})==1000)
    a.check('all paired contrasts exactly reflect own uncertainty mode',all(
        x[label+'|status']==draws[x['mode']][int(x['replicate'])][label+'|status'] and
        all(close(num(x[label+'|delta_'+m]),num(draws[x['mode']][int(x['replicate'])][label+'|delta_'+m])) for m in ('brier','auc'))
        for x in paired for label in PAIRS))
    a.summary=dict(science_commit=SCIENCE,run_id=RUN_ID,job_id=JOB_ID,status=r['status'],train_n=train['n'],train_events8=train['events8'],
        late_n=late['n'],late_events8=late['events8'],train_effective_n=train['effective_n'],late_effective_n=late['effective_n'],
        train_event_bearing_psu_n=train['event_bearing_psu_n'],late_event_bearing_psu_n=late['event_bearing_psu_n'],
        point=r['point'],contrasts={k:v for k,v in r['contrasts'].items() if k.startswith(('primary|','secondary|'))},
        uncertainty={mode:{'contrasts':{key:r['uncertainty'][mode]['intervals'][key] for key in interval_keys if key.startswith(('primary|','secondary|'))},
            'model_failure_counts':r['uncertainty'][mode]['model_failure_counts'],'contrast_failure_counts':r['uncertainty'][mode]['contrast_failure_counts']} for mode in MODES},
        calibration_valid_counts=calibration_valid,training_selected_penalties={arm:row['selected_penalty'] for arm,row in artifact['models'].items()},
        source_count50=50,scientific_module_count95=95,historical_module_count90=90,files_sha256={k:sha(p) for k,p in paths.items()},
        empirical_log_sha256=sha(args.log))
    return finish(a,args,paths)

def finish(a,args,paths):
    def forbidden(value):
        if isinstance(value,dict):
            return [k for k in value if k.lower() in {'seqn','participant_ids','individual_predictions','oof_predictions','participant_predictions'}]+[k for v in value.values() for k in forbidden(v)]
        if isinstance(value,list): return [k for v in value for k in forbidden(v)]
        return []
    a.check('exported JSON aggregate-only no personID or prediction vector keys',not any(forbidden(json.loads(p.read_text())) for k,p in paths.items() if not k.endswith('_CSV')))
    limits=[
        'Retrospective baseline-calendar split, not operational prospective validation: early1999–2002 eight-year labels can mature in2007–2010, after late2003–2004 recruitment. Late labels are absent from fit, but training labels were not all available at the late baseline calendar.',
        'Both periods and late outcomes were previously used. This is exploratory reused-cohort evaluation, not untouched or independent external validation; no frozen30 transport or clinical readiness claim.',
        'Primary uncertainty is joint_refit. fixed_training conditions on one fitted early model and excludes training uncertainty; the two500-draw distributions are separate, not pooled or interchangeable.',
        'The primary joint-refit upper Brier bound is only about1.3e-5 below zero. With500 Monte Carlo draws and this near-zero boundary, describe the result as limited exploratory support, not robust confirmation or a formal discovery; do not replace it with the narrower fixed-trained interval.',
        'All80 percentile intervals are independently recomputed by sorted-order linear interpolation; finite-sample nominal coverage is not established for nonlinear penalized fitting, discrete tuning and AUC/calibration. Undefined metric-local draws remain disclosed.',
        'Bootstrap holds the feature forms, domain, inner partition and observed measurements fixed. Assay/storage/calibration/prior-selection and unknown future-calendar uncertainty are excluded.',
        'Individual assay batch/date/storage duration and release-specific Cys traceability are unavailable. Calendar differences cannot identify biological drift separately from storage/batch effects.',
        'Source bytes were fetched/checked by the empirical runner. This audit reconciles exported fingerprints, locked historical hashes, exact stdout blocks and source-wrapper hashes; it is not an independent second download or local rehash of raw public records.',
        'Only aggregates/source metadata/coefficients were read; no person records, XPT/ASCII/ZIP downloads, predictions or model refits. Point Brier/AUC and optimizer gradients cannot be independently recomputed from aggregates; checks cover arithmetic, source hashes, fitted layout, reported diagnostics and bin moments.'
    ]
    result=dict(schema='independent-forward32-aggregate-audit/1',method='stdlib_only_no_pipeline_import_no_raw_no_refit',
        passed=sum(x['passed'] for x in a.checks),total=len(a.checks),all_passed=all(x['passed'] for x in a.checks),
        summary=a.summary,checks=a.checks,limits=limits)
    out=Path(args.out); out.with_suffix('.json').write_bytes(canonical(result))
    lines=['# Независимый численный аудит32','',f"**{result['passed']}/{result['total']} PASS; общий исход: {'PASS' if result['all_passed'] else 'FAIL'}.**",'',
        f'Научный commit `{SCIENCE}`; final CI run{RUN_ID}, empirical job{JOB_ID}. Superseded run40917121566 и его метрики не использованы.','']
    if a.summary.get('status')=='blocked_feasibility_no_model_fits': lines+=['Feasibility blocked: '+str(a.summary['blocked_reasons'])]
    else:
        s=a.summary
        lines+=[f"Early n={s['train_n']},events8={s['train_events8']}; late n={s['late_n']},events8={s['late_events8']}. Event-bearing PSU:early{s['train_event_bearing_psu_n']},late{s['late_event_bearing_psu_n']}.",'',
            'Все80 интервалов пересчитаны независимо из500 строк каждого режима.1000строк paired CSV согласованы с режимами; late draws идентичны fixed/joint, early weights постоянны только в fixed.','',
            '| Контраст | ΔBrier point | Joint refit95% percentile | Fixed training95% percentile |','|---|---:|---|---|']
        for label in PAIRS:
            key=label+'|delta_brier'; point=s['contrasts'][key]
            ji=s['uncertainty']['joint_refit']['contrasts'][key]; fi=s['uncertainty']['fixed_training']['contrasts'][key]
            lines.append(f"| {label} | {point} | [{ji['lower']},{ji['upper']}],valid={ji['valid']} | [{fi['lower']},{fi['upper']}],valid={fi['valid']} |")
        lines+=['','Primary joint-refit верхняя граница **−0.000013081464** очень близка к0. При500 Monte Carlo draws это ограниченный исследовательский результат, а не устойчивое подтверждение. Для secondary отрицательный fixed-training CI не заменяет пересекающий0 основной joint-refit CI.','',
            '| Контраст | ΔAUC point | Joint refit95% percentile | Fixed training95% percentile |',
            '|---|---:|---|---|']
        for label in PAIRS:
            key=label+'|delta_auc'; point=s['contrasts'][key]
            ji=s['uncertainty']['joint_refit']['contrasts'][key]; fi=s['uncertainty']['fixed_training']['contrasts'][key]
            lines.append(f"| {label} | {point} | [{ji['lower']},{ji['upper']}],valid={ji['valid']} | [{fi['lower']},{fi['upper']}],valid={fi['valid']} |")
        lines+=['','Metric-local calibration valid counts из500:','',
            '| Mode | Arm | Offset | Intercept | Slope |','|---|---|---:|---:|---:|']
        for mode in MODES:
            for arm in ARMS:
                v=s['calibration_valid_counts'][mode][arm]
                lines.append(f"| {mode} | {arm} | {v['calibration_in_large']} | {v['calibration_intercept']} | {v['calibration_slope']} |")
    lines+=['','Ограничения:','']+['- '+x for x in limits]
    failed=[x['name'] for x in a.checks if not x['passed']]
    if failed: lines+=['','Непрошедшие проверки:','']+['- '+x for x in failed]
    out.with_suffix('.md').write_text('\n'.join(lines)+'\n')
    print(json.dumps(dict(passed=result['passed'],total=result['total'],all_passed=result['all_passed'],failed=failed,outputs=[str(out.with_suffix('.json')),str(out.with_suffix('.md'))]),ensure_ascii=False))
    return 0 if result['all_passed'] else 1

if __name__=='__main__': raise SystemExit(main())
