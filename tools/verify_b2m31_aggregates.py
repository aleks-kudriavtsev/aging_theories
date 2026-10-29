#!/usr/bin/env python3
"""Independent aggregate-only audit of b2m31; stdlib, no model/data imports.

Prepared before empirical outputs. Does not fetch anything or fit any model.
Default evidence is the saved successful-job plaintext and repository exports.
"""
import argparse
import csv
import hashlib
import io
import json
import math
from pathlib import Path
import re
import subprocess
from collections import Counter
from urllib.parse import urlparse

ARMS = ('routine', 'routine_b2m', 'routine_cys', 'routine_cys_b2m')
PAIRS = {'primary': ('routine', 'routine_b2m'), 'secondary': ('routine_cys', 'routine_cys_b2m')}
METRICS = ('auc', 'brier', 'observed', 'predicted', 'bias_pp', 'observed_expected',
           'calibration_in_large', 'calibration_intercept', 'calibration_slope')
CORE = ['age_linear','age_quadratic','male','current','former','sbp','bmi','hba1c','egfr_low','egfr_high','log_uacr']
NAMES = ['age_centered','age_centered_squared','male','smoking_current','smoking_former','systolic_bp','bmi','hba1c','egfr_low','egfr_high','log2_uacr']
INDICES = {'routine': [], 'routine_b2m': [1], 'routine_cys': [0], 'routine_cys_b2m': [0,1]}
PNAMES = ['log2_cystatin_c','log2_beta2_microglobulin']
SCIENCE_COMMIT = 'f27c7b47157db7fbaf31dc998bbb5c55c1c4f7e8'
RUN_ID, JOB_ID = 40908831254, 123570909752

def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def finite(v):
    return isinstance(v, (int,float)) and not isinstance(v,bool) and math.isfinite(v)

def number(v):
    if v is None or v == '':
        return None
    n=float(v)
    if not math.isfinite(n):
        raise ValueError('Nonfinite exported numeric cell')
    return n

def close(a,b):
    return (a is None and b is None) or (finite(a) and finite(b) and math.isclose(a,b,rel_tol=1e-10,abs_tol=2e-12))

def percentile(values,p):
    """Independent order-statistic linear interpolation; no pipeline helper."""
    if not values:
        return None
    ordered=sorted(values)
    pos=(len(ordered)-1)*p
    lo=math.floor(pos); hi=math.ceil(pos)
    return ordered[lo]+(ordered[hi]-ordered[lo])*(pos-lo)

def canonical(obj):
    return (json.dumps(obj,ensure_ascii=False,indent=2,allow_nan=False)+'\n').encode()

def stripped_log(path):
    lines=Path(path).read_text().splitlines()
    return '\n'.join(re.sub(r'^\d{4}-\d\d-\d\dT\S+\s','',x) for x in lines)+'\n'

def extract(log,label):
    start='B2M31_'+label+'_BEGIN\n'; end='B2M31_'+label+'_END'
    if log.count(start)!=1 or log.count(end)!=1:
        raise ValueError('Expected exactly one complete marked block '+label)
    return log.split(start,1)[1].split(end,1)[0]

class Audit:
    def __init__(self):
        self.checks=[]
        self.summary={}
    def check(self,name,ok,detail=None):
        self.checks.append(dict(name=name,passed=bool(ok),**({'detail':detail} if detail is not None else {})))
    def group(self,name,conditions,detail=None):
        self.check(name,all(conditions),detail)

def point_check(a,point,n,events,wtotal,label):
    a.check(label+' arms',set(point)==set(ARMS))
    for arm in ARMS:
        p=point[arm]
        a.group(label+' '+arm+' population', [p['n']==n,p['events']==events,
            p['n_positive_weight']==n,p['events_positive_weight']==events,close(p['weight_sum'],wtotal)])
        a.group(label+' '+arm+' metric algebra',[
            all(finite(p[k]) and 0<=p[k]<=1 for k in ('auc','brier','observed','predicted')),
            close(p['bias_pp'],100*(p['predicted']-p['observed'])),
            close(p['observed_expected'],p['observed']/p['predicted'])])
        a.check(label+' '+arm+' calibration diagnostics',
            p['calibration_status']!='computed' or (all(finite(p[k]) for k in METRICS[-3:])
                and finite(p['calibration_gradient_max']) and p['calibration_gradient_max']<=1e-5))
    a.check(label+' common observed risk',all(close(point[x]['observed'],point[ARMS[0]]['observed']) for x in ARMS))

def contrast_check(a,point,contrast,label):
    for name,(base,expanded) in PAIRS.items():
        a.check(label+' '+name+' contrast status',contrast[name+'|status']=='computed')
        for metric in ('brier','auc'):
            a.check(label+' '+name+' delta '+metric,
                    close(contrast[name+'|delta_'+metric],point[expanded][metric]-point[base][metric]))
    a.check(label+' flattened arm values',all(close(contrast[arm+'|'+metric],point[arm][metric])
                                               for arm in ARMS for metric in METRICS))

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo',default=str(Path(__file__).resolve().parents[1]))
    parser.add_argument('--log',default='/tmp/b2m31_success_job.log')
    parser.add_argument('--ledger',default=str(Path(__file__).resolve().parents[1]/'research/mortality/b2m31/INDEPENDENT_EXPECTED_SOURCES.json'))
    parser.add_argument('--out',default='/tmp/verify31')
    parser.add_argument('--ci-jobs',default='/tmp/b2m31_independent_ci_jobs.json')
    parser.add_argument('--ci-evidence-dir',default='/tmp')
    args=parser.parse_args();ci_evidence=Path(args.ci_evidence_dir)
    repo=Path(args.repo); mortality=repo/'research/mortality'; root=mortality/'b2m31'
    a=Audit(); ledger=json.loads(Path(args.ledger).read_text()); lock=json.loads((root/'LOCK.json').read_text())
    spec=json.loads((root/'SOURCE_SPEC.json').read_text())
    log=stripped_log(args.log)
    paths={'RESULTS':root/'RESULTS_SUMMARY.json','MANIFEST':root/'SOURCE_MANIFEST.json',
           'PARAMETERS':root/'FOLD_PARAMETERS.json','CALIBRATION_BINS':root/'CALIBRATION_BINS.json',
           'BOOTSTRAP_CSV':root/'bootstrap_aggregates.csv','PRIMARY_BOOTSTRAP_CSV':root/'primary_bootstrap.csv',
           'SECONDARY_BOOTSTRAP_CSV':root/'secondary_bootstrap.csv'}
    if not paths['RESULTS'].exists():
        paths['RESULTS']=root/'RESULTS.json'
    if not paths['PARAMETERS'].exists():
        paths['PARAMETERS']=root/'fold_parameters.json'
    if not paths['CALIBRATION_BINS'].exists():
        paths['CALIBRATION_BINS']=root/'calibration_bins.json'
    data={}
    for label,p in paths.items():
        block=extract(log,label)
        if label.endswith('_CSV'):
            a.check(label+' exact bytes match successful log',p.read_bytes()==block.encode())
            data[label]=list(csv.DictReader(io.StringIO(block)))
        else:
            data[label]=json.loads(block)
            a.check(label+' saved structure matches successful log',json.loads(p.read_text())==data[label])
            a.check(label+' saved canonical bytes match successful log',p.read_bytes()==block.encode())
    r=data['RESULTS']; manifest=data['MANIFEST']; params=data['PARAMETERS']; bins=data['CALIBRATION_BINS']; sa=r['source_audit']
    wide=data['BOOTSTRAP_CSV']
    a.check('all500 refits completed in successful log','B2M31_NESTED_REFITS 500' in log)
    ci_jobs=json.loads(Path(args.ci_jobs).read_text())
    a.check('independently queried completed science CI jobs',len(ci_jobs)==3 and
        {j['name'] for j in ci_jobs}=={'empirical','synthetic (3.11)','synthetic (3.13)'} and
        all(j['run_id']==RUN_ID and j['status']=='completed' and j['conclusion']=='success' for j in ci_jobs))
    for j in ci_jobs:
        ownlog=stripped_log(ci_evidence/('b2m31_independent_job_'+str(j['id'])+'.log'))
        a.check(j['name']+' correct science commit included in CI checkout',
            'Merge '+SCIENCE_COMMIT+' into 50682455e6207785103dfc86fc7f57fd7b502a1f' in ownlog)
        if j['name'].startswith('synthetic'):
            a.check(j['name']+' 66 synthetic tests succeeded',bool(re.search(r'Ran 66 tests in [\d.]+s\n\s*\nOK\n',ownlog)))
        else:
            a.check('all seven blocks match independently fetched successful empirical log',
                all(extract(ownlog,label)==extract(log,label) for label in paths))
    a.check('analysis fixed before current metrics',r['plan_commit']=='6774fbe' and
            subprocess.check_output(['git','merge-base','--is-ancestor','6774fbe',SCIENCE_COMMIT],cwd=repo).decode()=='')
    a.group('prespecified fixed settings',[r['horizon_years']==8,r['horizon_months']==96,
       r['bootstrap_replicates']==500,r['bootstrap_seed']==20261102,
       r['primary_contrast']=='routine_b2m_minus_routine',r['secondary_contrast']=='routine_cys_b2m_minus_routine_cys',
       r['routine_features']==CORE,r['protein_features']==['SSCYST','SSB2M'],
       r['primary_score']=='survey_weighted_eight_year_Brier'])
    a.group('research scope explicit',[r['analysis']=='exploratory_reused_NHANES1999_2004_new_B2M_binary_models',
       r['independent_external_validation'] is False,r['clinical_use_ready'] is False,
       r['participant_records_exported'] is False,r['previous_frozen_coefficients_used'] is False,
       r['assay_calibration_and_prior_selection_uncertainty_included'] is False,
       r['bootstrap_partition_fixed'] is True,r['full_nested_preprocessing_tuning_and_refitting'] is True,
       r['same_common_domain_all_four_models'] is True,r['no_single_protein_winner_search'] is True,
       r['model_family']=='eight_year_survey_weighted_ridge_logistic_not_cause_specific_engine'])
    # All provenance checks use aggregate metadata or already-versioned modules.
    a.check('plan and lock hashes current',r['plan_sha256']==sha(root/'ANALYSIS_PLAN.md') and r['lock_sha256']==sha(root/'LOCK.json'))
    for rel,digest in lock['metadata_sha256'].items():
        p=root/rel
        if not p.exists(): p=mortality/rel
        a.check('metadata pin '+rel,p.exists() and sha(p)==digest)
    historical=lock['historical_scientific_python_sha256']
    a.check('84 historical scientific modules fixed',len(historical)==84 and all(sha(mortality/p)==v for p,v in historical.items()))
    current={str(p.relative_to(mortality)):sha(p) for p in sorted(mortality.rglob('*.py'))}
    a.check('complete90 current scientific module hashes',len(current)==90 and r['scientific_python_sha256']==current)
    a.check('six additions exclusively b2m31',len(set(current)-set(historical))==6 and all(p.startswith('b2m31/') for p in set(current)-set(historical)))
    a.check('manifest exact export SHA',sha(paths['MANIFEST'])==r['source_manifest_export_sha256'])
    a.check('manifest lock/schema/plan',[manifest.get(k) for k in ('schema','plan_commit','lock_sha256')]==
            ['b2m31-retrieval/1','6774fbe',sha(root/'LOCK.json')])
    expected={x['file']:x for x in ledger['required_file_inventory']}
    files={x['file']:x for x in manifest['files']}
    a.check('exact50 required sources unique',len(manifest['files'])==len(files)==50 and set(files)==set(expected))
    safe=lambda url: urlparse(url).scheme=='https' and urlparse(url).hostname in {'wwwn.cdc.gov','www.cdc.gov','ftp.cdc.gov'}
    a.check('required sources successful safe and nonempty',all(x['status']=='retrieved' and isinstance(x['bytes'],int) and x['bytes']>0
        and re.fullmatch('[0-9a-f]{64}',x['sha256']) and safe(x['url']) and safe(x['resolved_url']) for x in files.values()))
    a.check('fixed URLs match independent inventory',all(e['url'] is None or files[name]['url']==e['url'] for name,e in expected.items()))
    a.check('SSCARD dynamic official href filename',urlparse(files['1999_SSCARD_A.xpt']['url']).path.lower().endswith('/sscard_a.xpt'))
    a.check('25 historical public data SHA preserved',len(ledger['required_historical_data_sha256'])==25 and
            all(files[p]['sha256']==v for p,v in ledger['required_historical_data_sha256'].items()))
    denied=ledger['manifest_contract']['denied_documentation_outside_files']
    a.check('two403 PDF denials separately retained without retry',manifest['documentation_not_retrieved']==denied and sa['documentation_not_retrieved']==denied
            and all(x['file'] not in files and not any(k in x for k in ('sha256','bytes','resolved_url')) for x in denied))
    fingerprints={p:x['sha256'] for p,x in files.items()}
    a.check('source audit inventory and stable digest',sa['files_verified']==50 and sa['historical_data_hashes_verified']==25 and
            sa['source_fingerprints']==fingerprints and sa['source_content_sha256']==hashlib.sha256(json.dumps(fingerprints,sort_keys=True,separators=(',',':')).encode()).hexdigest())
    a.check('original runtime source manifest exact exported hash',sa['source_manifest_sha256']==sha(paths['MANIFEST']))
    # Full-release marginal expectations are codebook guards, not newly observed local records.
    expected_ss=ledger['full_SSCARD_release_checks']
    a.check('SSCARD all source marginal guards remain pinned',all(spec['sscard'][k]==expected_ss[k] for k in
        ('metadata_marginal_flag_counts','metadata_marginal_pristine_counts','metadata_weight_counts','expected_rows_full_release'))
        and sa['corrected_sscard_source_rows']==21890)
    a.check('corrected source weight rule unaltered',r['weight']=='source_qualified_six_year_weight' and
            '(2/3)*frame.WTSSCB4Y' in (root/'cohort.py').read_text() and '(1/3)*frame.WTSSCB2Y' in (root/'cohort.py').read_text())
    a.check('source-specific serum calibration',sa['serum_calibration_applied_by_cycle']=={'1999':'1.013*raw+0.147','2001':'no_additional_adjustment','2003':'no_additional_adjustment'})
    # Independent reader comparison hashes must agree with the already-pinned public release.
    readers=sa['mortality_reader_audits']
    refs={x['file']:x for x in ledger['mortality_reader_contract']['per_cycle_historical_reader_evidence']}
    a.check('three mortality readers exact8fields and249008cells',len(readers)==3 and
        sum(x['rows'] for x in readers)==31126 and sum(x['cell_comparisons'] for x in readers)==249008 and
        all(all(x[k]==refs[x['file']][k] for k in ('rows','fields','cell_comparisons','exact_match','sha256')) for x in readers))
    zero=sa['xpt_zero_audits']; zfiles={x['filename'] for x in zero}
    a.check('22 XPT zero-representation audits',len(zero)==len(zfiles)==22 and zfiles=={p for p in files if p.endswith('.xpt')}
        and all(x['sha256']==files[x['filename']]['sha256'] and x['epsilon_threshold_used'] is False
            and all(isinstance(v,int) and 0<=v<=x['nrows'] for v in x['numeric_zero_counts'].values())
            and all(0<=v<=x['numeric_zero_counts'][k] for k,v in x['corrected_counts'].items()) for x in zero))
    design=sa['full_DEMO_masked_design']; eligible_design=sa['eligible_masked_design']
    a.check('actual fullDEMO inventory exactly preserved at eligibility',design==eligible_design and sa['full_DEMO_design_inventory_preserved'] is True
        and sum(map(len,design.values()))==sa['full_DEMO_masked_psu_n']==sa['eligible_masked_psu_n'])
    distribution=dict(sorted(Counter(str(len(v)) for v in design.values()).items()))
    a.check('actual PSU inventory distribution and uniqueness',distribution==sa['eligible_psu_count_distribution'] and
            all(len(v)>=2 and len(v)==len(set(v)) for v in design.values()) and len(design)==43)
    a.check('PSU values match each cycle dictionary',all(any(lo<=int(st)<=hi and set(values)<=set(c['allowed_SDMVPSU'])
        for c in ledger['cycles'] for lo,hi in [c['SDMVSTRA_inclusive_range']]) for st,values in design.items()))
    flow=sa['cohort_flow']; common=sa['common_domain']; reported=sa['reported_sensitivity_domain']
    a.check('eligible flow reconciles cycles',len(flow)==3 and sum(x['eligible_weighted_subsample_n'] for x in flow)==r['eligible_n']==sa['eligible_n'] and
        all(x['survey_n']==next(c['full_DEMO_rows'] for c in ledger['cycles'] if c['year']==x['cycle']) and
            x['survey_n']>=x['age40_75_n']>=x['MEC_target_n']>=x['linked_target_n']>=x['eligible_weighted_subsample_n'] and
            x['linked_target_n']-x['eligible_weighted_subsample_n']==x['excluded_no_corrected_subsample_weight_n'] for x in flow))
    a.check('common domain n events effective n coverage',common['n']==r['n'] and common['events8']==r['events8'] and
        close(common['effective_n'],r['effective_n']) and close(sa['common_weighted_coverage'],r['common_weighted_coverage']) and
        0<r['effective_n']<=r['n']<=reported['n']<=r['eligible_n'] and 0<r['common_weighted_coverage']<=1)
    for label,domain,flags in [('primary',common,{'0.0'}),('reported sensitivity',reported,{'0.0','1.0'})]:
        a.group(label+' source domain marginals',[sum(domain['cycles'].values())==domain['n'],
            sum(domain['pristine_counts'].values())==domain['n'],0<domain['effective_n']<=domain['n'],
            all(sum(counts.values())==domain['n'] and set(counts)<=flags for counts in domain['flag_counts'].values()),
            all(domain['cycles'].get(str(x['cycle']),0)<=x['eligible_weighted_subsample_n'] for x in flow)])
    patterns=sa['exclusion_reason_patterns']; reasons=sa['exclusion_reason_counts']
    a.check('overlapping exclusions reconcile union and every reason',sum(patterns.values())==r['eligible_n'] and patterns.get('complete',0)==r['n'] and
        all(sum(n for p,n in patterns.items() if reason in p.split('|'))==v for reason,v in reasons.items()))
    a.check('routine missingness consistent with joint exclusion',set(sa['missing_routine'])==set(CORE) and
        max(sa['missing_routine'].values())<=reasons['routine']<=sum(sa['missing_routine'].values()))
    a.check('eightyear outcome and reuse limitations explicit',sa['minimum_survivor_followup_months']>=96 and
        sa['zero_month_deaths']>=0 and sa['source_all_participants_from_previously_loaded_release'] is True and
        sa['old_fitted_M4_participant_intersection_not_reconstructed'] is True and sa['prior_Cys_fitted'] is True and sa['prior_B2M_loaded'] is True)
    point_check(a,r['point'],r['n'],r['events8'],common['weight_sum'],'primary domain')
    contrast_check(a,r['point'],r['contrasts'],'primary domain')
    # Independent arithmetic over all500 wide aggregate rows, including undefined metric-local draws.
    a.check('all500 bootstrap identifiers ordered unique',[int(x['replicate']) for x in wide]==list(range(500)))
    interval_keys={arm+'|'+metric for arm in ARMS for metric in METRICS}|{label+'|delta_'+metric for label in PAIRS for metric in ('brier','auc')}
    a.check('exact40 prespecified interval targets',set(r['intervals'])==interval_keys)
    for key in sorted(interval_keys):
        vals=[number(x[key]) for x in wide if number(x[key]) is not None]
        interval=r['intervals'][key]
        a.group('independent percentile '+key,[close(interval['lower'],percentile(vals,.025)),close(interval['upper'],percentile(vals,.975)),
            interval['valid']==len(vals),interval['undefined']==500-len(vals),interval['defined_draws_only']==(len(vals)!=500),
            interval['scope']=='nested_preprocessing_tuning_refitting_with_fixed_partition'])
    wide_math=True
    for row in wide:
        for arm in ARMS:
            status=row[arm+'|status']; vals={m:number(row[arm+'|'+m]) for m in METRICS}
            wide_math &= status in {'computed','undefined'}
            if status=='undefined':
                wide_math &= bool(row[arm+'|error']) and all(v is None for v in vals.values())
            else:
                wide_math &= all(vals[k] is not None and 0<=vals[k]<=1 for k in ('brier','observed','predicted'))
                wide_math &= close(vals['bias_pp'],100*(vals['predicted']-vals['observed']))
                wide_math &= close(vals['observed_expected'],vals['observed']/vals['predicted'])
                if row[arm+'|calibration_status']=='computed':
                    wide_math &= all(vals[k] is not None for k in METRICS[-3:])
        for label,(base,expanded) in PAIRS.items():
            both=all(row[arm+'|status']=='computed' for arm in (base,expanded))
            wide_math &= row[label+'|status']==('computed' if both else 'undefined')
            for metric in ('brier','auc'):
                x,y=number(row[expanded+'|'+metric]),number(row[base+'|'+metric])
                wide_math &= close(number(row[label+'|delta_'+metric]),None if x is None or y is None else x-y)
    a.check('all500 arm metric/status and paired delta arithmetic',wide_math)
    for arm in ARMS:
        failed=sum(x[arm+'|status']!='computed' for x in wide)
        a.check(arm+' failure count and2500 fold selection counts',r['bootstrap_model_failure_counts'][arm]==failed and
            set(r['bootstrap_penalty_selection_counts'][arm])=={'1.0','10.0','100.0'} and
            sum(r['bootstrap_penalty_selection_counts'][arm].values())==5*(500-failed))
        a.check(arm+' calibration validity distinguished from model failure',all(
            r['intervals'][arm+'|'+metric]['valid']==sum(number(x[arm+'|'+metric]) is not None for x in wide)
            for metric in METRICS[-3:]) and
            sum(x[arm+'|calibration_status']=='computed' for x in wide)==
              r['intervals'][arm+'|calibration_slope']['valid'])
    for label in PAIRS:
        narrow=data[label.upper()+'_BOOTSTRAP_CSV']
        a.check(label+'500 narrow rows agree with full wide CSV',len(narrow)==500 and all(
            int(n['replicate'])==int(w['replicate']) and n['status']==w[label+'|status'] and
            all(close(number(n['delta_'+m]),number(w[label+'|delta_'+m])) for m in ('brier','auc')) for n,w in zip(narrow,wide)))
        a.check(label+' paired contrast failure count',r['bootstrap_contrast_failure_counts'][label]==sum(w[label+'|status']!='computed' for w in wide))
    for label,key in [('BOOTSTRAP_CSV','bootstrap_aggregates_csv_sha256'),('PRIMARY_BOOTSTRAP_CSV','primary_bootstrap_csv_sha256'),('SECONDARY_BOOTSTRAP_CSV','secondary_bootstrap_csv_sha256')]:
        a.check(label+' exact file hash',sha(paths[label])==r[key])
    sensitivity=r['reported_values_sensitivity']
    a.check('sensitivity distinct conditional domain no bootstrap pooling',sensitivity['n']==reported['n'] and
        sensitivity['interpretation']=='reported_values_flags01_different_conditional_domain_within_domain_pairs_only' and sensitivity['bootstrap_intervals_computed'] is False)
    domains={'primary_common':(r['n'],r['events8'],r['point'],common['weight_sum'])}
    if sensitivity['status']=='computed':
        point_check(a,sensitivity['point'],reported['n'],reported['events8'],reported['weight_sum'],'reported sensitivity')
        contrast_check(a,sensitivity['point'],sensitivity['contrasts'],'reported sensitivity')
        a.check('sensitivity event count',sensitivity['events8']==reported['events8'])
        domains['reported_values_sensitivity']=(reported['n'],reported['events8'],sensitivity['point'],reported['weight_sum'])
    else:
        a.check('undefined sensitivity fit honestly retained',sensitivity['status']=='undefined' and bool(sensitivity.get('error')))
    a.check('exact5fold fourarm parameter records per computed domain',len(params)==20*len(domains) and
        Counter((p['analysis_domain'],p['model']) for p in params)==Counter({(d,arm):5 for d in domains for arm in ARMS}))
    for domain,(n,events,point,wtotal) in domains.items():
        rows=[p for p in params if p['analysis_domain']==domain]
        for arm in ARMS:
            fs=[p for p in rows if p['model']==arm]; dim=11+len(INDICES[arm])
            a.check(domain+' '+arm+' fold totals',sorted(p['outer_fold'] for p in fs)==list(range(5)) and
                sum(p['scoring_n'] for p in fs)==n and sum(p['scoring_events'] for p in fs)==events and
                all(p['training_n']+p['scoring_n']==n and p['training_events']+p['scoring_events']==events and
                    p['training_positive_weight_n']==p['training_n'] and p['training_positive_weight_events']==p['training_events'] and
                    p['PSU_overlap_n']==0 and p['training_psu_n']>=3 and p['scoring_psu_n']>=1 for p in fs))
            a.check(domain+' '+arm+' complete feature/scaling layout',all(p['core_feature_count']==11 and p['core_feature_indices']==list(range(11)) and
                p['protein_feature_indices']==INDICES[arm] and p['feature_names']==NAMES+[PNAMES[j] for j in INDICES[arm]] and
                len(p['core_coefficients_standardized'])==11 and len(p['protein_coefficients_standardized'])==len(INDICES[arm]) and
                len(p['training_mean'])==len(p['training_scale'])==len(p['constant_features'])==dim and
                all(finite(v) for v in p['training_mean']+p['training_scale']+p['core_coefficients_standardized']+p['protein_coefficients_standardized']) and
                all(v>0 for v in p['training_scale']) and all(isinstance(v,bool) for v in p['constant_features']) and
                all(not const or scale==1. for const,scale in zip(p['constant_features'],p['training_scale'])) and
                p['protein_transformation']=='reported_log2' for p in fs))
            a.check(domain+' '+arm+' independent penalty selection',all(set(p['inner_brier_by_penalty'])=={'1.0','10.0','100.0'} and
                all(finite(v) and 0<=v<=1 for v in p['inner_brier_by_penalty'].values()) and
                p['lambda']==p['selected_penalty']==max(float(k) for k,v in p['inner_brier_by_penalty'].items() if v<=min(p['inner_brier_by_penalty'].values())+1e-12) for p in fs))
            a.check(domain+' '+arm+' convergence/normalization',all(finite(p['intercept']) and finite(p['objective']) and p['iterations']>0 and
                close(p['normalized_weight_sum'],p['training_n']) and p['training_weight_sum']>0 and
                finite(p['gradient_max']) and ((p['convergence_status']=='gradient_1e-6' and p['gradient_max']<=1e-6) or
                    (p['convergence_status']=='tiny_step_gradient_1e-4' and p['gradient_max']<=1e-4)) for p in fs))
            bs=[b for b in bins if b['analysis_domain']==domain and b['model']==arm]
            a.check(domain+' '+arm+' calibration-bin counts',1<=len(bs)<=5 and len({b['bin'] for b in bs})==len(bs) and
                sum(b['n'] for b in bs)==n and sum(b['events'] for b in bs)==events and
                sum(b['n_positive_weight'] for b in bs)==n and close(sum(b['weight_sum'] for b in bs),wtotal))
            a.check(domain+' '+arm+' independent calibration-bin weighted moments',all(b['weight_sum']>0 and 0<=b['observed']<=1 and
                0<=b['predicted']<=1 and close(b['bias_pp'],100*(b['predicted']-b['observed'])) for b in bs) and
                close(sum(b['weight_sum']*b['observed'] for b in bs)/wtotal,point[arm]['observed']) and
                close(sum(b['weight_sum']*b['predicted'] for b in bs)/wtotal,point[arm]['predicted']))
        foldkeys=('training_n','scoring_n','training_events','scoring_events','training_psu_n','scoring_psu_n','training_weight_sum')
        a.check(domain+' same wholePSU folds across all four arms',all(len({tuple(p[k] for k in foldkeys) for p in rows if p['outer_fold']==fold})==1 for fold in range(5)))
        a.check(domain+' training/scoring totalweights reconcile',all(close(sum(p['training_weight_sum'] for p in rows if p['model']==arm),4*wtotal) for arm in ARMS))
        if domain=='primary_common':
            a.check('primary total PSU count from each fold',all(p['training_psu_n']+p['scoring_psu_n']==r['psu_n'] for p in rows))
    a.check('calibration bins exclusively computed domains/arms',all(b['analysis_domain'] in domains and b['model'] in ARMS for b in bins))
    # Explicitly reject accidental identifier or individual prediction vector keys.
    def badkeys(value):
        if isinstance(value,dict):
            return [k for k,v in value.items() if k.lower() in {'seqn','participant_ids','participant_predictions','oof_predictions','individual_predictions'}]+[k for v in value.values() for k in badkeys(v)]
        if isinstance(value,list): return [k for v in value for k in badkeys(v)]
        return []
    a.check('aggregate exports have no individual ID or prediction vector keys',not any(badkeys(data[k]) for k in ('RESULTS','MANIFEST','PARAMETERS','CALIBRATION_BINS')))
    a.summary=dict(science_commit=SCIENCE_COMMIT,run_id=RUN_ID,job_id=JOB_ID,eligible_n=r['eligible_n'],n=r['n'],events8=r['events8'],
        effective_n=r['effective_n'],strata_n=r['strata_n'],psu_n=r['psu_n'],eligible_psu_n=sa['eligible_masked_psu_n'],
        eligible_psu_distribution=sa['eligible_psu_count_distribution'],primary_delta_brier=r['contrasts']['primary|delta_brier'],
        primary_delta_auc=r['contrasts']['primary|delta_auc'],secondary_delta_brier=r['contrasts']['secondary|delta_brier'],
        secondary_delta_auc=r['contrasts']['secondary|delta_auc'],primary_brier_interval=r['intervals']['primary|delta_brier'],
        secondary_brier_interval=r['intervals']['secondary|delta_brier'],model_failures=r['bootstrap_model_failure_counts'],
        paired_failures=r['bootstrap_contrast_failure_counts'],sensitivity_status=sensitivity['status'],sensitivity_n=sensitivity['n'],
        folds=len(params),scientific_modules=len(current),historical_modules=len(historical),sources=len(files),
        historical_data_sources=25,mortality_reader_cell_comparisons=249008,files_sha256={k:sha(v) for k,v in paths.items()},
        successful_job_log_sha256=sha(args.log),point=r['point'],
        calibration_valid_counts={arm:{metric:r['intervals'][arm+'|'+metric]['valid'] for metric in METRICS[-3:]} for arm in ARMS},
        calibration_status_counts={arm:dict(Counter(x[arm+'|calibration_status'] for x in wide)) for arm in ARMS},
        synthetic_tests_per_python=66,independently_queried_ci_jobs=ci_jobs,
        ci_log_sha256={j['name']:sha(ci_evidence/('b2m31_independent_job_'+str(j['id'])+'.log')) for j in ci_jobs})
    limits=[
      'No raw participant/XPT/ASCII/ZIP downloads, individual prediction access or refitting occurred in this audit.',
      'Source SHA values are reconciled across locked historical pins, current exported manifest, runtime audit and success log; this is not an independent local rehash of raw source bytes.',
      'Full-release SSCARD marginal counts and corrected weight counts are checked by the pinned runtime guard; no separate observed full-release marginal table is exported. The independent audit verifies those guards and exported conditional-domain marginals, not a second XPT parse.',
      'Person-level Brier/AUC and fitted calibration regressions cannot be recomputed from aggregate exports; independently checked quantities are interval arithmetic, pair differences, metric identities and calibration-bin weighted moments.',
      '500 rescaled-stratified-PSU bootstrap draws repeat training, scaling and inner tuning but hold the original whole-PSU partitions fixed. Assay standardization and prior model/marker selection uncertainty are excluded.',
      'Cohort1999-2004 and Cys have prior use; B2M was previously loaded. These are exploratory new binary models on reused persons, not independent external validation.',
      'Both official mortality PDFs remained denied HTTP403, denied_not_retried, separate from50 required sources. The amendment retains prior primary definition review, pinned public ASCII hashes and two-reader agreement without changing the endpoint.',
      'The flags0/1 sensitivity uses a different conditional domain and supplies within-domain paired point comparisons only; it is not an additional bootstrap interval or external validation.'
    ]
    result=dict(schema='independent-b2m31-aggregate-audit/1',method='stdlib_only_no_pipeline_import_no_raw_no_refit',
        passed=sum(x['passed'] for x in a.checks),total=len(a.checks),all_passed=all(x['passed'] for x in a.checks),
        summary=a.summary,checks=a.checks,limits=limits)
    out=Path(args.out)
    out.with_suffix('.json').write_bytes(canonical(result))
    failed=[x['name'] for x in a.checks if not x['passed']]
    report=['# Независимый численный аудит итерации31','',f"Проверки: **{result['passed']}/{result['total']}**; общий исход: **{'PASS' if result['all_passed'] else 'FAIL'}**.",'',
       f'Научный commit: `{SCIENCE_COMMIT}`; CI run{RUN_ID}, job{JOB_ID}. Независимый запрос GitHub подтвердил empirical success и66 synthetic tests на каждой Python3.11/3.13. CI checkout — merge ref с этим commit; научные хеши совпадают с опубликованными результатами.',
       f"Общий домен: n={r['n']}, событий8лет={r['events8']}, effective n={r['effective_n']:.6f}; eligible n={r['eligible_n']}; {r['strata_n']}strata/{r['psu_n']}PSU.",'',
       'Аудитор использует только стандартную библиотеку и сохранённые агрегаты. Проверены точные marked-blocks, в том числе по независимо полученному completed-job log,50 источников и25 исторических data SHA,12 metadata pins,84 сохранённых и6 новых научных модулей,500 строк каждого CSV, все40 интервалов, парная арифметика, параметры фолдов и калибровочные корзины.','',
       '| Контраст | ΔBrier | percentile95% CI | ΔAUC | percentile95% CI |','|---|---:|---|---:|---|']
    for label in PAIRS:
        db=r['contrasts'][label+'|delta_brier']; da=r['contrasts'][label+'|delta_auc']; ib=r['intervals'][label+'|delta_brier']; ia=r['intervals'][label+'|delta_auc']
        report.append(f"| {label} | {db:.10f} | [{ib['lower']}, {ib['upper']}] valid={ib['valid']} | {da:.10f} | [{ia['lower']}, {ia['upper']}] valid={ia['valid']} |")
    report+=['',f"Fold records: {len(params)}; sensitivity: {sensitivity['status']}, n={sensitivity['n']}. Model failures: `{r['bootstrap_model_failure_counts']}`.",'',
       'У калибровки меньше определённых bootstrap-реплик, чем у Brier/AUC:','',
       '| Модель | Offset valid/500 | Intercept valid/500 | Slope valid/500 |',
       '|---|---:|---:|---:|']
    for arm in ARMS:
        v=a.summary['calibration_valid_counts'][arm]
        report.append(f"| {arm} | {v['calibration_in_large']} | {v['calibration_intercept']} | {v['calibration_slope']} |")
    report+=['','Эти калибровочные интервалы условны по определённым диагностическим fits. Нулевые ошибки прогнозных refits не означают500 определённых slope/intercept.','',
       'Ограничения:','']+['- '+x for x in limits]
    if failed: report+=['','Непрошедшие проверки:','']+['- '+x for x in failed]
    out.with_suffix('.md').write_text('\n'.join(report)+'\n')
    print(json.dumps(dict(passed=result['passed'],total=result['total'],all_passed=result['all_passed'],failed=failed,outputs=[str(out.with_suffix('.json')),str(out.with_suffix('.md'))]),ensure_ascii=False))
    return 0 if result['all_passed'] else 1

if __name__=='__main__':
    raise SystemExit(main())
