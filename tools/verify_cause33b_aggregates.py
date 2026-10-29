#!/usr/bin/env python3
"""Independent stdlib audit helpers for cause33B, before predictive outputs.

No scientific-pipeline imports, raw parsing, network requests, fitting or repo
writes. Final schema adapter will be pinned after the protocol is committed.
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

PLAN_COMMIT='d772b874198a9f7739225bbce5d1b6b6dfcd0c08'
BOOTSTRAP_SEED=20261105

CAUSES=('002','001')
ARMS=('routine','routine_b2m','routine_cys','routine_cys_b2m')
METRICS=('auc','brier','observed','predicted','bias_pp','observed_expected',
    'calibration_in_large','calibration_intercept','calibration_slope')
PAIRS={'primary':('routine','routine_b2m'),'secondary':('routine_cys','routine_cys_b2m')}
SOURCE_GATE_RESULTS_SHA='13c8d9b42f510ffee6e7b04476e545b02f9496ec08f53b585263950badf180b5'
SOURCE_GATE_MANIFEST_SHA='a52c16a4bf7bf48b581301aa638b0978849d1122454b0ca2bff3fdbb1abf88e1'
CORE=['age_linear','age_quadratic','male','current','former','sbp','bmi','hba1c','egfr_low','egfr_high','log_uacr']
CORE_NAMES=['age_centered','age_centered_squared','male','smoking_current','smoking_former','systolic_bp','bmi','hba1c','egfr_low','egfr_high','log2_uacr']
PNAMES=['log2_cystatin_c','log2_beta2_microglobulin']
INDICES={'routine':[],'routine_b2m':[1],'routine_cys':[0],'routine_cys_b2m':[0,1]}
SCIENCE_FIELDS={'schema','gate','endpoint','eligible','common','per_cycle','causes','unknown_cause',
    'defined_cause_ranking','priority_top5','proceed_codes','accounting','design','folds','weight_rule',
    'fitted_models','predictive_metrics_computed','participant_records_exported','threshold_scope'}

def finite(value):
    return isinstance(value,(int,float)) and not isinstance(value,bool) and math.isfinite(value)

def close(a,b):
    return (a is None and b is None) or (finite(a) and finite(b) and math.isclose(a,b,rel_tol=1e-10,abs_tol=2e-12))

def number(value):
    if value is None or value=='':return None
    result=float(value)
    if not math.isfinite(result):raise ValueError('Nonfinite aggregate numeric value')
    return result

def percentile(values,p):
    """Independent order-statistic interpolation, no numpy/pipeline helper."""
    if not values:return None
    ordered=sorted(values);position=(len(ordered)-1)*p
    low=math.floor(position);high=math.ceil(position)
    return ordered[low]+(ordered[high]-ordered[low])*(position-low)

def summarize_interval(values,total=500):
    valid=[v for v in values if v is not None and finite(v)]
    return dict(lower=percentile(valid,.025),upper=percentile(valid,.975),
        valid=len(valid),undefined=total-len(valid),defined_draws_only=len(valid)!=total)

def paired_delta(baseline,expanded):
    return None if baseline is None or expanded is None else expanded-baseline

def hashed_folds(keys,salt,folds):
    ordered=sorted(set(keys),key=lambda key:hashlib.sha256((salt+':'+key).encode()).hexdigest())
    return {key:i%folds for i,key in enumerate(ordered)}

def map_digest(value):
    return hashlib.sha256(json.dumps(value,sort_keys=True,separators=(',',':')).encode()).hexdigest()

def registered_interval_keys(cause_prefix=True):
    local=[arm+'|'+metric for arm in ARMS for metric in METRICS]
    local += [label+'|delta_'+metric for label in PAIRS for metric in ('brier','auc')]
    return [cause+'|'+key for cause in CAUSES for key in local] if cause_prefix else local

def source_control_counts(n,allcause_events,cause_events,unknown_early=0):
    if unknown_early:
        raise ValueError('Early unknown causes block cause33B comparisons')
    if not 0<=cause_events<=allcause_events<=n:
        raise ValueError('Invalid disjoint source counts')
    return dict(event_rows=cause_events,known_noncase_rows=n-cause_events,
        known_competing_early_deaths=allcause_events-cause_events,
        no_death_by_horizon_rows=n-allcause_events)

def metric_algebra(point):
    return close(point['bias_pp'],100*(point['predicted']-point['observed'])) and \
        close(point['observed_expected'],point['observed']/point['predicted'] if point['predicted']>0 else None)

def scientific_projection_compare(actual,reference):
    """Independent strict JSON/count comparator with only fixed float tolerance."""
    result=dict(floating_values_checked=0,floating_values_different=0,maximum_float_abs_difference=0.)
    def visit(a,b):
        if type(a) is not type(b):raise ValueError('Scientific projection type mismatch')
        if isinstance(b,dict):
            if set(a)!=set(b):raise ValueError('Scientific projection field mismatch')
            for key in b:visit(a[key],b[key])
        elif isinstance(b,list):
            if len(a)!=len(b):raise ValueError('Scientific projection length mismatch')
            for aa,bb in zip(a,b):visit(aa,bb)
        elif type(b) is float:
            if not math.isfinite(a) or not math.isfinite(b) or not math.isclose(a,b,rel_tol=1e-12,abs_tol=1e-9):
                raise ValueError('Scientific floating reduction outside locked tolerance')
            result['floating_values_checked']+=1
            result['floating_values_different']+=int(a!=b)
            result['maximum_float_abs_difference']=max(result['maximum_float_abs_difference'],abs(a-b))
        elif type(b) in (str,int,bool,type(None)):
            if a!=b:raise ValueError('Scientific count/label/fold mismatch')
        else:raise ValueError('Scientific projection contains non-JSON type')
    visit(actual,reference)
    return result

def canonical(value):
    return (json.dumps(value,ensure_ascii=False,indent=2,allow_nan=False)+'\n').encode()

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def log_text(path):
    # GitHub's joined log chunks may carry BOM before an ISO timestamp even
    # mid-log. Remove it only together with that exact transport prefix.
    prefix=r'^\ufeff?\d{4}-\d\d-\d\dT\d\d:\d\d:\d\d(?:\.\d+)?Z ?'
    return '\n'.join(re.sub(prefix,'',line) for line in Path(path).read_text().splitlines())+'\n'

def marked(log,name):
    start='CA33B_FILE_BEGIN '+name+'\n';end='CA33B_FILE_END '+name
    if log.count(start)!=1 or log.count(end)!=1:raise ValueError('Nonunique/incomplete33B marked file '+name)
    return log.split(start,1)[1].split(end,1)[0]

class Audit:
    def __init__(self):self.checks=[];self.summary={}
    def check(self,name,passed,detail=None):
        self.checks.append(dict(name=name,passed=bool(passed),**({'detail':detail} if detail is not None else {})))

def converged(record):
    g=record['gradient_max']
    return finite(g) and ((record['convergence_status']=='gradient_1e-6' and g<=1e-6) or
        (record['convergence_status']=='tiny_step_gradient_1e-4' and g<=1e-4))

def metric_row_math(row):
    """Local arm/pair/metric validity, without asserting AUC defined on fit success."""
    ok=True
    for arm in ARMS:
        status=row[arm+'|status'];values={m:number(row[arm+'|'+m]) for m in METRICS}
        ok &= status in {'computed','undefined'}
        if status=='undefined':
            ok &= bool(row.get(arm+'|error')) and all(v is None for v in values.values())
        else:
            ok &= all(values[k] is not None and 0<=values[k]<=1 for k in ('brier','observed','predicted'))
            ok &= values['auc'] is None or 0<=values['auc']<=1
            ok &= metric_algebra(values)
            if row[arm+'|calibration_status']=='computed':ok &= all(values[k] is not None for k in METRICS[-3:])
    for label,(base,expanded) in PAIRS.items():
        both=all(row[arm+'|status']=='computed' for arm in (base,expanded))
        ok &= row[label+'|status']==('computed' if both else 'undefined')
        for metric in ('brier','auc'):
            ok &= close(number(row[label+'|delta_'+metric]),paired_delta(number(row[base+'|'+metric]),number(row[expanded+'|'+metric])))
    return ok

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo',default='.')
    parser.add_argument('--log',default='/tmp/cause33B_success_job.log')
    parser.add_argument('--out',default='/tmp/cause33B_audit')
    parser.add_argument('--science-commit',required=True)
    parser.add_argument('--run-id',required=True,type=int)
    parser.add_argument('--job-id',required=True,type=int)
    parser.add_argument('--tests-per-python',required=True,type=int)
    parser.add_argument('--ci-jobs',default='/tmp/cause33B_independent_ci_jobs.json')
    parser.add_argument('--ci-log-dir',default='/tmp')
    parser.add_argument('--ci-log-prefix',default='cause33B_independent_job_')
    args=parser.parse_args();repo=Path(args.repo);mortality=repo/'research/mortality';root=mortality/'cause33b'
    a=Audit();log=log_text(args.log);lock=json.loads((root/'LOCK.json').read_text());contract=json.loads((root/'CONTRACT.json').read_text())
    names=['RESULTS_SUMMARY.json','FOLD_PARAMETERS.json','CALIBRATION_BINS.json','SOURCE_MANIFEST.json','bootstrap_aggregates.csv','paired_bootstrap.csv']
    files={name:root/name for name in names};data={}
    for name,path in files.items():
        text=marked(log,name)
        a.check(name+' exact completed-job marked bytes',path.read_bytes()==text.encode())
        data[name]=list(csv.DictReader(io.StringIO(text))) if name.endswith('.csv') else json.loads(text)
    r=data['RESULTS_SUMMARY.json'];params=data['FOLD_PARAMETERS.json'];bins=data['CALIBRATION_BINS.json'];manifest=data['SOURCE_MANIFEST.json']
    wide=data['bootstrap_aggregates.csv'];paired=data['paired_bootstrap.csv']
    a.check('protocol committed before final sciencehead',subprocess.run(['git','merge-base','--is-ancestor',PLAN_COMMIT,args.science_commit],cwd=repo,capture_output=True).returncode==0)
    a.check('plan contract and lock hashes',r['plan_commit']==PLAN_COMMIT and r['plan_sha256']==sha(root/'ANALYSIS_PLAN.md') and
        r['contract_sha256']==sha(root/'CONTRACT.json') and r['lock_sha256']==sha(root/'LOCK.json'))
    for name,digest in lock['metadata_sha256'].items():
        path=root/name if name in lock['local_metadata'] else mortality/name
        a.check('metadata pin '+name,sha(path)==digest)
    actual={str(p.relative_to(mortality)):sha(p) for p in sorted(mortality.rglob('*.py'))};old=lock['historical_scientific_python_sha256']
    a.check('all99 historical scientific modules unchanged',len(old)==99 and all(actual.get(k)==v for k,v in old.items()))
    a.check('all103 scientific module hashes current',len(actual)==103 and r['scientific_python_sha256']==actual)
    a.check('four new modules restricted to cause33b',len(set(actual)-set(old))==4 and all(k.startswith('cause33b/') for k in set(actual)-set(old)))
    for name,digest in r['outputs_sha256'].items():a.check('exact result output SHA '+name,sha(files[name])==digest)
    a.check('source manifest export SHA',r['source_manifest_export_sha256']==sha(files['SOURCE_MANIFEST.json']))
    a.check('exact fixed cases controls common domain horizon',[r['causes'],r['n'],r['eligible_n'],r['horizon_months']]==[list(CAUSES),6044,6546,96] and
        r['events8']=={'002':176,'001':148} and r['controls8']=={'002':5868,'001':5896} and
        r['routine_features']==CORE and r['protein_features']==['SSCYST','SSB2M'])
    a.check('fixed bootstrap and models from contract',r['bootstrap_replicates']==contract['bootstrap_replicates']==500 and
        r['bootstrap_seed']==contract['bootstrap_seed']==BOOTSTRAP_SEED and contract['causes']==list(CAUSES) and contract['arms']==list(ARMS) and
        contract['contrasts']=={k:list(v) for k,v in PAIRS.items()} and contract['common_n']==6044 and contract['allcause_events8']==586)
    a.check('correct controls multiplicity and scope flags',r['auc_controls']==contract['auc_controls']=='all_known_noncases_including_early_competing_deaths_late_deaths_survivors_observed_through96' and
        r['weight']=='pooled_corrected_six_year_SSCARD_weight' and r['marginal_exploratory_intervals_not_familywise_confirmation'] is True and
        r['between_cause_effect_comparisons_performed'] is False and r['independent_binary_probabilities_summed'] is False and
        r['joint_integration_and_Fine_Gray']=='deferred_not_implemented_or_run' and r['independent_external_validation'] is False and
        r['clinical_use_ready'] is False and r['participant_records_exported'] is False and
        r['original_models_and_lab_report_v3_presentation_changed'] is False and r['comparative_model_fitting_attempted'] is True)
    a.check('paired fullrefit fixed domain partition declarations',r['shared_preprocessing_verified_for_available_matched_fits'] is True and
        r['shared_preprocessing_verification_scope']=='available_same_arm_outer_fold_pairs_only; unavailable_pairs_after_fit_failure_not_verified' and
        r['one_eligible_PSU_draw_per_replicate_shared_across_causes_and_arms'] is True and r['retained_zero_weight_rows_and_fixed_folds'] is True and r['bootstrap_redraw'] is False)
    a.check('all500 attempted full-refits completed','CA33B_SHARED_FULL_REFITS 500' in log)
    # Independently compare the scientific source projection to the frozen33A gate.
    sa=r['source_audit'];inherited=sa['inherited31_source_domain_audit'];reproduced=sa['reproduced33A_gate'];proof=sa['scientific_gate_equality']
    frozen=json.loads((mortality/'cause33/GATE_RESULTS.json').read_text())
    reference={key:frozen[key] for key in SCIENCE_FIELDS};reference['schema']='cause33-count-gate/1'
    a.check('frozen33A source pins exact',sha(mortality/'cause33/GATE_RESULTS.json')==SOURCE_GATE_RESULTS_SHA==sa['frozen33A_gate_sha256']==manifest['frozen33A_gate_sha256'] and
        sha(mortality/'cause33/SOURCE_MANIFEST.json')==SOURCE_GATE_MANIFEST_SHA)
    projection_stats=scientific_projection_compare(reproduced,reference)
    a.check('source33A science projection exact types counts fixed tolerance',set(reproduced)==SCIENCE_FIELDS and proof['exact_structure_and_counts'] is True and
        proof['numeric_tolerance']=={'rtol':1e-12,'atol':1e-9} and proof['expected_proceed_codes']==list(CAUSES) and
        proof['compared_fields']==sorted(SCIENCE_FIELDS) and all(proof[k]==v for k,v in projection_stats.items()) and
        proof['recomputed_scientific_sha256']==map_digest(reproduced) and proof['frozen_scientific_sha256']==map_digest(reference))
    a.check('source stage proves two supported targets before fitting',reproduced['proceed_codes']==list(CAUSES) and
        reproduced['gate']['source_endpoint_qualified'] is True and reproduced['common']['unknown_cause_deaths8']==0 and
        all(reproduced['causes'][cause]['proceed_to_later_protocol'] is True for cause in CAUSES) and
        sa['comparative_models_fitted'] is False and sa['labels_exported'] is False and sa['participant_records_exported'] is False)
    base=manifest['source'];intake=manifest['intake'];fingerprints={x['file']:x['sha256'] for x in base['files']}
    a.check('source wrapper current plan and base hashes',manifest['schema']=='cause33b-source-manifest/1' and intake['schema']=='cause33b-retrieval/1' and
        intake['plan_commit']==PLAN_COMMIT and intake['lock_sha256']==sha(root/'LOCK.json') and
        intake['base_manifest_sha256']==hashlib.sha256(canonical(base)).hexdigest()==inherited['source_manifest_sha256'])
    a.check('all50 source fingerprints linked',len(base['files'])==len(fingerprints)==50 and intake['source_fingerprints']==fingerprints==inherited['source_fingerprints'] and
        inherited['files_verified']==50 and inherited['historical_data_hashes_verified']==25 and inherited['source_content_sha256']==map_digest(fingerprints))
    old_data={name:digest for digest,name in (line.split() for line in (mortality/'nhanes04/SOURCE_DATA_SHA256.txt').read_text().splitlines() if line.strip())}
    safe=lambda url:urlparse(url).scheme=='https' and urlparse(url).hostname in {'wwwn.cdc.gov','www.cdc.gov','ftp.cdc.gov'}
    a.check('official retrieval rows and25 historical data SHA',all(x['status']=='retrieved' and x['bytes']>0 and safe(x['url']) and safe(x['resolved_url']) for x in base['files']) and
        len([name for name in fingerprints if name.endswith(('.dat','.xpt'))])==25 and all(old_data.get(name)==digest for name,digest in fingerprints.items() if name.endswith(('.dat','.xpt'))))
    denied=base['documentation_not_retrieved']
    a.check('two previousPDF403 never retried or replaced',len(denied)==2 and denied==inherited['documentation_not_retrieved'] and all(
        x['status']=='denied_not_retried' and x['original_http_status']==403 and x['file'] not in fingerprints and 'sha256' not in x for x in denied))
    a.check('uncensoring and two-reader mortality retained',inherited['minimum_survivor_followup_months']>=96 and
        len(inherited['mortality_reader_audits'])==3 and sum(x['cell_comparisons'] for x in inherited['mortality_reader_audits'])==249008 and
        all(x['exact_match'] is True and x['fields']==8 and x['sha256']==fingerprints[x['file']] for x in inherited['mortality_reader_audits']))
    total_weight=reproduced['common']['weight_sum'];global_n=r['n']
    # All8 point models use the same full control set, including competing deaths.
    for cause in CAUSES:
        a.check(cause+' allcause gate disjoint control accounting',source_control_counts(global_n,586,r['events8'][cause])==dict(
            event_rows=r['events8'][cause],known_noncase_rows=r['controls8'][cause],
            known_competing_early_deaths=586-r['events8'][cause],no_death_by_horizon_rows=5458))
        a.check(cause+' point pair and local metric arithmetic',metric_row_math(r['contrasts'][cause]))
        for arm in ARMS:
            point=r['point'][cause].get(arm);status=r['point_statuses'][cause][arm]
            a.check(cause+' '+arm+' point status retained',r['contrasts'][cause][arm+'|status']==status['status'] and
                ((point is not None)==(status['status']=='computed')))
            if point is None:continue
            a.check(cause+' '+arm+' identical n control weight denominator',point['n']==global_n and point['events']==r['events8'][cause] and
                point['n_positive_weight']==global_n and point['events_positive_weight']==r['events8'][cause] and close(point['weight_sum'],total_weight) and
                close(point['observed'],reproduced['causes'][cause]['overall']['weighted_events8']/total_weight) and
                all(close(point[m],r['contrasts'][cause][arm+'|'+m]) for m in METRICS))
            a.check(cause+' '+arm+' calibration diagnostic convergence',point['calibration_status']!='computed' or (
                finite(point['calibration_gradient_max']) and point['calibration_gradient_max']<=1e-5 and all(finite(point[m]) for m in METRICS[-3:])))
            bs=[b for b in bins if b['cause']==cause and b['model']==arm]
            a.check(cause+' '+arm+' calibration-bin counts and weighted moments',1<=len(bs)<=5 and len({b['bin'] for b in bs})==len(bs) and
                sum(b['n'] for b in bs)==global_n and sum(b['events'] for b in bs)==r['events8'][cause] and sum(b['n_positive_weight'] for b in bs)==global_n and
                all(b['weight_sum']>0 and close(b['bias_pp'],100*(b['predicted']-b['observed'])) for b in bs) and close(sum(b['weight_sum'] for b in bs),total_weight) and
                close(sum(b['weight_sum']*b['observed'] for b in bs)/total_weight,point['observed']) and
                close(sum(b['weight_sum']*b['predicted'] for b in bs)/total_weight,point['predicted']))
    # Point training artifacts:40records if all8 predictive fits complete.
    expected_fits=sum(not (r['point_statuses'][cause][arm]['status']=='undefined' and r['point_statuses'][cause][arm].get('error','').startswith('fit|')) for cause in CAUSES for arm in ARMS)
    a.check('five pointfold records per completed nested fit',len(params)==expected_fits*5 and all(p['cause'] in CAUSES and p['model'] in ARMS for p in params))
    for cause in CAUSES:
        for arm in ARMS:
            fs=[p for p in params if p['cause']==cause and p['model']==arm]
            if not fs:continue
            e=r['events8'][cause];dim=11+len(INDICES[arm])
            a.check(cause+' '+arm+' pointfold totals no PSU leakage',sorted(p['outer_fold'] for p in fs)==list(range(5)) and sum(p['scoring_n'] for p in fs)==global_n and
                sum(p['scoring_events'] for p in fs)==e and all(p['training_n']+p['scoring_n']==global_n and p['training_events']+p['scoring_events']==e and
                    p['training_positive_weight_n']==p['training_n'] and p['training_positive_weight_events']==p['training_events'] and p['PSU_overlap_n']==0 and
                    p['training_psu_n']+p['scoring_psu_n']==87 for p in fs) and close(sum(p['training_weight_sum'] for p in fs),4*total_weight))
            a.check(cause+' '+arm+' registered features scales coefficients',all(p['core_feature_count']==11 and p['core_feature_indices']==list(range(11)) and
                p['protein_feature_indices']==INDICES[arm] and p['feature_names']==CORE_NAMES+[PNAMES[i] for i in INDICES[arm]] and
                len(p['training_mean'])==len(p['training_scale'])==len(p['constant_features'])==dim and
                len(p['core_coefficients_standardized'])==11 and len(p['protein_coefficients_standardized'])==len(INDICES[arm]) and
                all(finite(v) for v in p['training_mean']+p['training_scale']+p['core_coefficients_standardized']+p['protein_coefficients_standardized']) and
                all(v>0 for v in p['training_scale']) and all(not c or scale==1 for c,scale in zip(p['constant_features'],p['training_scale'])) and
                p['protein_transformation']=='reported_log2' for p in fs))
            a.check(cause+' '+arm+' independent innerBrier lambda tie selection',all(set(p['inner_brier_by_penalty'])=={'1.0','10.0','100.0'} and
                all(finite(v) and 0<=v<=1 for v in p['inner_brier_by_penalty'].values()) and p['lambda']==p['selected_penalty']==max(
                    float(k) for k,v in p['inner_brier_by_penalty'].items() if v<=min(p['inner_brier_by_penalty'].values())+1e-12) for p in fs))
            a.check(cause+' '+arm+' gradients normalizedtraining likelihood',all(finite(p['intercept']) and finite(p['objective']) and p['iterations']>0 and
                converged(p) and close(p['normalized_weight_sum'],p['training_n']) for p in fs))
    invariant=['training_mean','training_scale','constant_features','feature_names','core_feature_indices','protein_feature_indices',
        'training_n','scoring_n','training_positive_weight_n','training_weight_sum','training_psu_n','scoring_psu_n','PSU_overlap_n']
    for arm in ARMS:
        for fold in range(5):
            rows=[p for p in params if p['model']==arm and p['outer_fold']==fold]
            a.check(arm+' fold'+str(fold)+' exact crosscause preprocessing and frame identity',len(rows)<2 or all(rows[0][k]==rows[1][k] for k in invariant))
    # Independent fixed mapping matches source gate, point PSU counts and all arms.
    inventory=inherited['full_DEMO_masked_design'];keys=[st+':'+str(p) for st,psus in inventory.items() for p in psus]
    outer=hashed_folds(keys,'b2m31:outer',5)
    a.check('independent preserved31 outermap fingerprint',reproduced['folds']['outer_mapping_sha256']==map_digest(outer) and len(keys)==87)
    for row in reproduced['folds']['outer']:
        f=row['outer_fold'];inner=hashed_folds([k for k,v in outer.items() if v!=f],'b2m31:inner:'+str(f),3)
        a.check('outer'+str(f)+' independent innerfold fingerprint',row['inner_mapping_sha256']==map_digest(inner) and
            all(p['scoring_psu_n']==sum(v==f for v in outer.values()) for p in params if p['outer_fold']==f))
    # Shared500wide draws: evaluate80 intervals independently, never treat the2causes as independent.
    a.check('all500 ordered retained wide attempts',[int(row['replicate']) for row in wide]==list(range(500)))
    a.check('all2000 paired causecontrast tuples unique',len(paired)==2000 and len({(int(p['replicate']),p['cause'],p['contrast']) for p in paired})==2000 and
        Counter((p['cause'],p['contrast']) for p in paired)==Counter({(c,k):500 for c in CAUSES for k in PAIRS}))
    a.check('shared draw valid positive weights and zero-row retention',all(number(x['eligible_weight_sum'])>=number(x['common_weight_sum'])>0 and
        number(x['eligible_weight_squared_sum'])>=number(x['common_weight_squared_sum'])>0 and
        0<number(x['common_positive_weight_rows'])<=global_n and number(x['common_positive_weight_rows'])<=number(x['eligible_positive_weight_rows'])<=6546 and
        0<number(x['common_positive_weight_PSUs'])<=87 and number(x['common_weight_sum'])**2/number(x['common_weight_squared_sum'])<=number(x['common_positive_weight_rows'])+1e-8 and
        all(re.fullmatch('[0-9a-f]{64}',x[k]) for k in ('shared_common_weight_sha256','shared_eligible_weight_sha256')) for x in wide))
    a.check('500variable shared-draw weight fingerprints',len({x['shared_common_weight_sha256'] for x in wide})>1 and len({x['shared_eligible_weight_sha256'] for x in wide})>1)
    a.check('disjoint cause positive-event weights no larger than samecommon frame',all(sum(number(x[c+'|event_weight_sum']) for c in CAUSES)<=number(x['common_weight_sum'])+1e-8 and
        sum(number(x[c+'|positive_weight_event_rows']) for c in CAUSES)<=number(x['common_positive_weight_rows']) for x in wide))
    calibration_valid={}
    for cause in CAUSES:
        local=[{key[len(cause)+1:]:value for key,value in row.items() if key.startswith(cause+'|')} for row in wide]
        a.check(cause+'500 arm pair status and metric identities',all(metric_row_math(row) for row in local))
        a.check(cause+' same shared draw denominator and knownnoncase controls',all(
            number(w[cause+'|positive_weight_event_rows'])+number(w[cause+'|positive_weight_noncase_rows'])==number(w['common_positive_weight_rows']) and
            0<=number(w[cause+'|positive_weight_event_rows'])<=r['events8'][cause] and
            all(row[arm+'|status']!='computed' or close(number(row[arm+'|observed']),number(w[cause+'|event_weight_sum'])/number(w['common_weight_sum'])) for arm in ARMS)
            for row,w in zip(local,wide)))
        a.check(cause+' exactly40 registered interval targets',set(r['intervals'][cause])==set(registered_interval_keys(False)))
        for key in registered_interval_keys(False):
            vals=[number(row[key]) for row in local];expected=summarize_interval(vals);actual_i=r['intervals'][cause][key]
            a.check(cause+' independent percentile '+key,all(close(actual_i[k],v) if k in ('lower','upper') else actual_i[k]==v for k,v in expected.items()) and
                actual_i['scope']=='nested_preprocessing_tuning_refitting_with_fixed_partition' and actual_i['marginal_exploratory'] is True and
                actual_i['simultaneous_familywise_coverage'] is False and actual_i['endpoint_selection_uncertainty_included'] is False)
        failures=r['bootstrap_failures'][cause]
        for arm in ARMS:
            a.check(cause+' '+arm+' separate fit-metric calibration failure reasons',failures['model_failure_counts'][arm]==sum(row[arm+'|status']!='computed' for row in local) and
                failures['fit_or_metric_error_reasons'][arm]==dict(Counter(row[arm+'|error'] for row in local if row[arm+'|error'])) and
                failures['calibration_undefined_reasons'][arm]==dict(Counter(row[arm+'|calibration_status'] or 'predictive_fit_or_metric_unavailable' for row in local if row[arm+'|calibration_status']!='computed')))
            fit_fail=sum(row[arm+'|status']!='computed' and row[arm+'|error'].startswith('fit|') for row in local)
            selection=r['bootstrap_penalty_selection_counts'][cause][arm]
            a.check(cause+' '+arm+'5foldlambda counts persuccessful predictive fit',set(selection)=={'1.0','10.0','100.0'} and all(isinstance(v,int) and v>=0 for v in selection.values()) and sum(selection.values())==5*(500-fit_fail))
        for contrast in PAIRS:a.check(cause+' '+contrast+' paired failure count',failures['contrast_failure_counts'][contrast]==sum(row[contrast+'|status']!='computed' for row in local))
        calibration_valid[cause]={arm:{m:r['intervals'][cause][arm+'|'+m]['valid'] for m in METRICS[-3:]} for arm in ARMS}
    a.check('all2000 paired values and metric-local statuses exactly matchwide',all(p['cause'] in CAUSES and p['contrast'] in PAIRS and
        p['status']==wide[int(p['replicate'])][p['cause']+'|'+p['contrast']+'|status'] and
        all(close(number(p['delta_'+m]),number(wide[int(p['replicate'])][p['cause']+'|'+p['contrast']+'|delta_'+m])) and
            p[m+'_status']==('computed' if number(p['delta_'+m]) is not None else 'undefined') for m in ('brier','auc')) for p in paired))
    # CI evidence is saved separately from the completed final job, with no network inside this script.
    if Path(args.ci_jobs).exists():
        jobs=json.loads(Path(args.ci_jobs).read_text())
        a.check('independent final CI all3jobs success',len(jobs)==3 and all(j['run_id']==args.run_id and j['status']=='completed' and j['conclusion']=='success' for j in jobs))
        for job in jobs:
            ownlog=log_text(Path(args.ci_log_dir)/(args.ci_log_prefix+str(job['id'])+'.log'))
            a.check(job['name']+' scientific commit in checkout',args.science_commit in ownlog)
            if job['name'].startswith('synthetic'):
                a.check(job['name']+' exact synthetic testcount',bool(re.search(r'Ran '+str(args.tests_per_python)+r' tests in [\d.]+s\n\s*\nOK\n',ownlog)))
            else:a.check('six exports match independently fetched completed empiricallog',all(marked(ownlog,name)==marked(log,name) for name in names))
    else:a.check('completed independent CI evidence supplied',False)
    def forbidden(value):
        if isinstance(value,dict):return [k for k in value if k.lower() in {'seqn','participant_ids','individual_labels','individual_predictions','participant_predictions','oof_predictions','bootstrap_coefficients'}]+[k for v in value.values() for k in forbidden(v)]
        if isinstance(value,list):return [k for v in value for k in forbidden(v)]
        return []
    a.check('aggregate exports lack participant keys labels probabilities orbootstrapcoefficients',not any(forbidden(data[name]) for name in names if name.endswith('.json')))
    merge=re.search(r'git log -1 --format=%H\n([0-9a-f]{40})',log)
    a.summary=dict(science_commit=args.science_commit,plan_commit=PLAN_COMMIT,run_id=args.run_id,job_id=args.job_id,
        actual_PR_merge_checkout=merge.group(1) if merge else None,tests_per_python=args.tests_per_python,
        n=global_n,events8=r['events8'],controls8=r['controls8'],point=r['point'],
        contrasts={cause:{key:value for key,value in r['contrasts'][cause].items() if key.startswith(('primary|','secondary|'))} for cause in CAUSES},
        contrast_intervals={cause:{key:r['intervals'][cause][key] for key in registered_interval_keys(False) if key.startswith(('primary|','secondary|'))} for cause in CAUSES},
        calibration_valid_counts=calibration_valid,bootstrap_failures=r['bootstrap_failures'],point_fit_records=len(params),
        scientific_files=len(actual),historical_files=len(old),source_files=50,data_files=25,source_projection=projection_stats,
        output_sha256={name:sha(path) for name,path in files.items()},empirical_log_sha256=sha(args.log))
    limits=[
        'Exploratory, internally nested comparisons on reused persons/proteins after full-cohort outcome-based source/support selection; no independent external or prospective validation.',
        'All80 intervals are marginal conditional percentile summaries, without familywise coverage, multiplicity-adjusted confirmation, p-values or a between-cause interaction contrast. Different interval signs or sizes do not establish different marker effects between causes.',
        'Controls include known early competing deaths, deathsafter96months and observed-through96 survivors. A competing death is a noncase of the fixed-horizon target, not censoring or proof of being alive.',
        'One eligible-frame PSU draw is reused across causes/arms; weight SHA values and moment/observed-risk invariants are checked, but the unpublished weight vector cannot be independently reconstructed or rehashed from aggregates.',
        '500 refits repeat preprocessing/tuning/training with fixed rows/folds/forms. Metric-local undefined draws are retained; defined-only bounds are not complete nominal95% coverage. Assay/storage/calibration, public perturbation, endpoint/model-selection history and transport uncertainty are excluded.',
        'Source hashes, source projection, point parameter layout/diagnostics, bin moments, deltas and all interval percentiles are independently reconciled; no raw person-record read or model refit occurred. Person-level Brier/AUC and optimizer gradients are not recomputed from aggregates.',
        'Independent binary cause probabilities are not summed or advertised as a coherent joint model. Fine–Gray/joint integration, frozen models and laboratory reportv3/presentation remain outside this change.'
    ]
    report=dict(schema='independent-cause33B-aggregate-audit/1',passed=sum(x['passed'] for x in a.checks),total=len(a.checks),
        all_passed=all(x['passed'] for x in a.checks),summary=a.summary,checks=a.checks,limits=limits)
    out=Path(args.out);out.with_suffix('.json').write_bytes(canonical(report))
    lines=['# Независимый численный аудит33B','',f"**{report['passed']}/{report['total']} PASS; общий исход: {'PASS' if report['all_passed'] else 'FAIL'}.**",'',
        f"Protocol `{PLAN_COMMIT}`; scientific commit `{args.science_commit}`; CIrun{args.run_id},job{args.job_id}; actual merge `{a.summary['actual_PR_merge_checkout']}`.",
        f"Synthetic Python3.11/3.13: {args.tests_per_python} tests each. Domain n={global_n}; events002={r['events8']['002']},001={r['events8']['001']}; knowncontrols5868/5896.",'',
        'Проверены80marginalintervals,500общихwide rows,2000pairedrows,source33Aprojection,99старыхмодулей,40pointfoldrecords при8успешныхfits иcrosscausepreprocessing invariance.','',
        '| Cause | Contrast | ΔBrier | Marginal95% percentile | ΔAUC | Marginal95% percentile |','|---|---|---:|---|---:|---|']
    for cause in CAUSES:
        for label in PAIRS:
            db=r['contrasts'][cause][label+'|delta_brier'];da=r['contrasts'][cause][label+'|delta_auc'];ib=r['intervals'][cause][label+'|delta_brier'];ia=r['intervals'][cause][label+'|delta_auc']
            lines.append(f"| {cause} | {label} | {db} | [{ib['lower']},{ib['upper']}], valid={ib['valid']} | {da} | [{ia['lower']},{ia['upper']}], valid={ia['valid']} |")
    lines+=['','Это не familywise confirmation и не interaction test; разный знак/определённость по причинам не доказывает различия эффектов.','',
        '| Cause | Arm | Offset valid/500 | Intercept valid/500 | Slope valid/500 |','|---|---|---:|---:|---:|']
    for cause in CAUSES:
        for arm in ARMS:
            v=calibration_valid[cause][arm]
            lines.append(f"| {cause} | {arm} | {v['calibration_in_large']} | {v['calibration_intercept']} | {v['calibration_slope']} |")
    lines+=['','Ограничения:','']+['- '+x for x in limits]
    failed=[x['name'] for x in a.checks if not x['passed']]
    if failed:lines+=['','Failed checks:','']+['- '+x for x in failed]
    out.with_suffix('.md').write_text('\n'.join(lines)+'\n')
    print(json.dumps(dict(passed=report['passed'],total=report['total'],all_passed=report['all_passed'],failed=failed),ensure_ascii=False))
    return 0 if report['all_passed'] else 1

if __name__=='__main__':raise SystemExit(main())
