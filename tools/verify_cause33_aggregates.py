#!/usr/bin/env python3
"""Independent stdlib cause33A audit helpers, prepared before empirical counts.

No pipeline imports, network calls, raw parsing, model fitting or repo writes.
Final GATE_RESULTS schema adapter is added after gate helper is finalized.
"""
import argparse
from collections import Counter
import hashlib
import json
import math
from pathlib import Path
import re
import subprocess
from urllib.parse import urlparse

PLAN_COMMIT='82eacda12eced3c0105019fee09236b647821f45'
DEFINED_CODES=tuple(f'{i:03d}' for i in range(1,10))
ALL_CODES=tuple(f'{i:03d}' for i in range(1,11))
THRESHOLDS={
 'event_rows_min':50,'noncase_rows_min':50,'event_psus_min':20,
 'outer_evaluation_event_rows_min':5,'outer_evaluation_noncase_rows_min':5,
 'outer_training_event_rows_min':20,'outer_training_noncase_rows_min':20,
 'inner_evaluation_event_rows_min':3,'inner_evaluation_noncase_rows_min':3,
 'inner_training_event_rows_min':10,'inner_training_noncase_rows_min':10}

def finite(value):
    return isinstance(value,(int,float)) and not isinstance(value,bool) and math.isfinite(value)

def close(a,b):
    return (a is None and b is None) or (finite(a) and finite(b) and math.isclose(a,b,rel_tol=1e-10,abs_tol=2e-12))

def hashed_folds(keys,salt,folds):
    """Independent SHA256 ordering; public masked-PSU keys, never person IDs."""
    ordered=sorted(set(keys),key=lambda key:hashlib.sha256((salt+':'+key).encode()).hexdigest())
    return {key:rank%folds for rank,key in enumerate(ordered)}

def rank_defined(records):
    """Canonical records map code -> weighted_events8/raw_events8."""
    return sorted(DEFINED_CODES,key=lambda code:(-records[code]['weighted_events8'],-records[code]['raw_events8'],code))

def burden_arithmetic(event_weight,total_weight,allcause_event_weight,squared_weight,max_event_weight):
    return dict(weighted_incidence8=event_weight/total_weight if total_weight>0 else None,
        weighted_death_share8=event_weight/allcause_event_weight if allcause_event_weight>0 else None,
        event_weight_kish_n=event_weight**2/squared_weight if squared_weight>0 else None,
        max_single_event_weight_fraction=max_event_weight/event_weight if event_weight>0 else None)

def support_failures(global_counts,cells):
    """Independent canonical failure list; adapter maps final JSON fields only.

    Cell kind: outer_training/outer_evaluation/inner_training/inner_evaluation.
    Counts always distinguish unknown early causes from known noncases.
    """
    failures=[]
    for key,rule in [('event_rows','event_rows_min'),('known_noncase_rows','noncase_rows_min'),('event_psus','event_psus_min')]:
        if global_counts[key]<THRESHOLDS[rule]:
            failures.append(dict(context='global',rule=rule,actual=global_counts[key],minimum=THRESHOLDS[rule]))
    for cell in cells:
        context=cell['kind']+':outer'+str(cell['outer'])
        if cell['kind'].startswith('inner_'):context+=':inner'+str(cell['inner'])
        for key,suffix in [('event_rows','event_rows_min'),('known_noncase_rows','noncase_rows_min')]:
            rule=cell['kind']+'_'+suffix
            if cell[key]<THRESHOLDS[rule]:
                failures.append(dict(context=context,rule=rule,actual=cell[key],minimum=THRESHOLDS[rule]))
    return failures

def known_noncase(n,event_rows,unknown_early_rows):
    result=n-event_rows-unknown_early_rows
    if min(n,event_rows,unknown_early_rows,result)<0:
        raise ValueError('Impossible disjoint endpoint accounting')
    return result

def canonical(value):
    return (json.dumps(value,ensure_ascii=False,indent=2,allow_nan=False)+'\n').encode()

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def digest(value):
    return hashlib.sha256(value).hexdigest()

def map_digest(value):
    return digest(json.dumps(value,sort_keys=True,separators=(',',':')).encode())

def log_text(path):
    return '\n'.join(re.sub(r'^\d{4}-\d\d-\d\dT\S+\s','',line)
        for line in Path(path).read_text().splitlines())+'\n'

def marked(log,filename):
    start='CAUSE33_STDOUT_FILE_BEGIN '+filename+'\n'; end='CAUSE33_STDOUT_FILE_END '+filename
    if log.count(start)!=1 or log.count(end)!=1:
        raise ValueError('Expected exactly one complete marked block '+filename)
    return log.split(start,1)[1].split(end,1)[0]

class Audit:
    def __init__(self):
        self.checks=[]; self.summary={}
    def check(self,name,passed,detail=None):
        self.checks.append(dict(name=name,passed=bool(passed),**({'detail':detail} if detail is not None else {})))

def cell_key(cell):
    return (cell['outer_fold'],cell['level'],cell['subset'],cell.get('inner_fold'))

def failure_key(failure):
    return (failure['rule'],failure['required'],failure['observed'],failure.get('outer_fold'),
        failure.get('level'),failure.get('subset'),failure.get('inner_fold'))

def independently_expected_failures(overall,cells):
    failures=[]
    for rule,value in [('event_rows_min',overall['positive_weight_event_rows']),
        ('noncase_rows_min',overall['noncase_rows']),('event_psus_min',overall['event_PSU_n'])]:
        if value<THRESHOLDS[rule]: failures.append(dict(rule=rule,required=THRESHOLDS[rule],observed=value))
    for cell in cells:
        context={key:cell[key] for key in ('outer_fold','inner_fold','level','subset') if key in cell}
        for suffix,value in [('event_rows_min',cell['positive_weight_event_rows']),('noncase_rows_min',cell['noncase_rows'])]:
            rule=cell['level']+'_'+cell['subset']+'_'+suffix
            if value<THRESHOLDS[rule]:failures.append(dict(rule=rule,required=THRESHOLDS[rule],observed=value,**context))
    return failures

def check_summary(a,label,value,population,clinical=True,unknown_n=None):
    e=value['events8']; ew=value['weighted_events8']; square=value['event_weight_squared_sum']; maximum=value['max_event_weight8']
    a.check(label+' event counts and label accounting',value['n']==population['n'] and
        value['positive_weight_event_rows']==e and 0<=e<=value['n'] and
        ((value['label_unknown_rows'] is None and value['noncase_rows'] is None) if not clinical else
         (value['label_unknown_rows']==unknown_n and value['noncase_rows']==known_noncase(value['n'],e,unknown_n))))
    a.check(label+' event weights Kish concentration',
        finite(ew) and ew>=0 and finite(square) and square>=0 and
        ((ew==square==0 and maximum is None and value['event_Kish_n'] is None and value['maximum_event_weight_fraction'] is None) if e==0 else
         (ew>0 and square>0 and finite(maximum) and 0<maximum<=ew and
          close(value['event_Kish_n'],ew**2/square) and 0<value['event_Kish_n']<=e+1e-9 and
          close(value['maximum_event_weight_fraction'],maximum/ew))))
    a.check(label+' survey burden arithmetic',close(value['common_weight_sum'],population['weight_sum']) and
        close(value['weighted_incidence8'],ew/population['weight_sum']) and
        close(value['weighted_share_of_all_deaths8'],ew/population['weighted_all_cause_events8'] if population['weighted_all_cause_events8'] else None) and
        close(value['raw_share_of_all_deaths8'],e/population['all_cause_events8'] if population['all_cause_events8'] else None))
    a.check(label+' event-bearing design bounds',0<=value['event_strata_n']<=value['event_PSU_n']<=min(87,e) and
        (e==0 or 0<value['event_strata_n']<=43))

def population_check(a,label,population):
    d=population
    a.check(label+' all-time death and source-cause accounting',d['alive_n']+d['all_time_deaths_n']==d['n'] and
        sum(d['all_time_cause_counts'].values())==d['all_time_deaths_n'] and
        set(d['all_time_cause_counts'])==set(ALL_CODES)|{'unknown'} and
        d['all_time_cause_counts']['unknown']==d['unknown_cause_deaths8']+d['unknown_cause_deaths_after96'] and
        0<=d['unknown_cause_deaths8']<=d['all_cause_events8']<=d['all_time_deaths_n']<=d['n'] and d['weight_sum']>0)

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--repo',default='.')
    p.add_argument('--log',default='/tmp/cause33A_success_job.log')
    p.add_argument('--out',default='/tmp/verify33A')
    p.add_argument('--science-commit',required=True)
    p.add_argument('--run-id',type=int,required=True)
    p.add_argument('--job-id',type=int,required=True)
    p.add_argument('--ci-jobs',default='/tmp/cause33A_independent_ci_jobs.json')
    p.add_argument('--ci-log-dir',default='/tmp')
    p.add_argument('--ci-log-prefix',default='cause33_independent_job_')
    p.add_argument('--tests-per-python',type=int,required=True)
    args=p.parse_args(); repo=Path(args.repo); mortality=repo/'research/mortality'; root=mortality/'cause33'
    a=Audit(); log=log_text(args.log)
    paths={name:root/name for name in ('GATE_RESULTS.json','SOURCE_MANIFEST.json')}
    data={}
    for name,path in paths.items():
        text=marked(log,name)
        a.check(name+' exact completed-job marked bytes',path.read_bytes()==text.encode())
        data[name]=json.loads(text)
    r=data['GATE_RESULTS.json']; manifest=data['SOURCE_MANIFEST.json']; lock=json.loads((root/'LOCK.json').read_text())
    rules=json.loads((root/'SUPPORT_RULES.json').read_text()); definitions=json.loads((root/'CAUSE_DEFINITIONS.json').read_text())
    a.check('completed count-only gate marker','CAUSE33_GATE_COMPLETE no_comparative_fits' in log)
    a.check('protocol precedes final scientific commit',subprocess.run(['git','merge-base','--is-ancestor',PLAN_COMMIT,args.science_commit],cwd=repo,capture_output=True).returncode==0)
    a.check('qualified schema and plan-lock binding',r['schema']=='cause33-qualified-counts/1' and r['plan_commit']==PLAN_COMMIT and r['lock_sha256']==sha(root/'LOCK.json'))
    for name,expected in lock['metadata_sha256'].items():
        path=root/name if (root/name).is_file() else mortality/name
        a.check('metadata pin '+name,sha(path)==expected)
    actual={str(path.relative_to(mortality)):sha(path) for path in sorted(mortality.rglob('*.py'))}
    old=lock['historical_scientific_python_sha256']
    a.check('95 historical scientific files unchanged',len(old)==r['historical_scientific_files_verified']==95 and all(actual.get(k)==v for k,v in old.items()))
    a.check('all current scientific hashes match gate runner',r['scientific_python_sha256']==actual)
    a.check('four new modules exclusively cause33',len(actual)==99 and len(set(actual)-set(old))==4 and all(k.startswith('cause33/') for k in set(actual)-set(old)))
    a.check('all exact precommitted support thresholds',all(rules[k]==v for k,v in THRESHOLDS.items()) and rules['rank_maximum']==5 and rules['rank_before_support_not_replace_failed_top5'] is True)
    a.check('fixed underlying-cause definition and no coarsened cycles',definitions['underlying_cause_field']=='UCOD_LEADING' and definitions['supported_baseline_start_years']==[1999,2001,2003] and
        definitions['horizon_months']==96 and set(definitions['groups'])==set(ALL_CODES) and all(
            definitions['groups'][code]['eligible_for_defined_clinical_ranking']==(code!='010') for code in ALL_CODES))
    a.check('stageA does not fit or claim marker effects',r['fitted_models'] is False and r['comparison_models_fitted'] is False and
        r['predictive_metrics_computed'] is False and r['participant_records_exported'] is False and
        r['independent_binary_probabilities_summed'] is False and r['joint_integration_and_Fine_Gray']=='deferred_not_implemented_or_run' and
        r['prior_cohort_use']=='all_cycles_previously_used_no_independent_validation')
    a.check('endpoint separates competing and unknown causes',r['endpoint']==dict(horizon_months=96,underlying_cause_field='UCOD_LEADING',
        competing_deaths_are_known_noncase=True,unknown_early_deaths_are_not_noncase=True))
    # Source provenance is inspected as aggregate metadata, without downloading source records.
    base=manifest['base_manifest']; wrapper=manifest['wrapper']; source=manifest['inherited_source_domain_audit']
    a.check('source wrapper binds exact base/protocol/lock',manifest['schema']=='cause33-qualified-sources/1' and wrapper['schema']=='cause33-retrieval/1' and
        wrapper['plan_commit']==PLAN_COMMIT and wrapper['lock_sha256']==sha(root/'LOCK.json') and
        wrapper['base_manifest_sha256']==digest(canonical(base))==source['source_manifest_sha256'])
    fingerprints={row['file']:row['sha256'] for row in base['files']}
    a.check('all50 official input entries and fingerprints',len(base['files'])==len(fingerprints)==50 and
        wrapper['source_fingerprints']==source['source_fingerprints']==fingerprints and source['files_verified']==50 and
        source['source_content_sha256']==map_digest(fingerprints))
    safe=lambda url:urlparse(url).scheme=='https' and urlparse(url).hostname in {'wwwn.cdc.gov','www.cdc.gov','ftp.cdc.gov'}
    a.check('50 retrieved nonempty official sources',all(row['status']=='retrieved' and row['bytes']>0 and safe(row['url']) and safe(row['resolved_url']) for row in base['files']))
    historical={name:hashvalue for hashvalue,name in (line.split() for line in (mortality/'nhanes04/SOURCE_DATA_SHA256.txt').read_text().splitlines() if line.strip())}
    a.check('25 old public data fingerprints preserved',source['historical_data_hashes_verified']==25 and len([k for k in fingerprints if k.endswith(('.dat','.xpt'))])==25 and
        all(historical.get(name)==value for name,value in fingerprints.items() if name.endswith(('.dat','.xpt'))))
    denied=base['documentation_not_retrieved']
    a.check('two priorPDF403 denials separate without retry or fictitious hashes',len(denied)==2 and denied==source['documentation_not_retrieved'] and all(
        x['status']=='denied_not_retried' and x['original_http_status']==403 and x['file'] not in fingerprints and 'sha256' not in x for x in denied))
    a.check('retained endpoint contract/review honest source-byte status',manifest['retained_endpoint_contract_sha256']==sha(mortality/'endpoint_contract11.py') and
        manifest['endpoint_review_sha256']==sha(root/'ENDPOINT_REVIEW.md') and manifest['dictionary_status']=='official_indexed_text_and_retained_review_not_new_PDF_bytes' and
        manifest['participant_rows_exported'] is False and manifest['ZIP_downloaded'] is False)
    a.check('current mortality readers agree3cycles249008cells',len(source['mortality_reader_audits'])==3 and
        sum(x['rows'] for x in source['mortality_reader_audits'])==31126 and sum(x['cell_comparisons'] for x in source['mortality_reader_audits'])==249008 and
        all(x['exact_match'] is True and x['fields']==8 and x['sha256']==fingerprints[x['file']] for x in source['mortality_reader_audits']))
    common=r['common']; eligible=r['eligible']; cycles=r['per_cycle']; causes=r['causes']; unknown=r['unknown_cause']
    a.check('exact inherited pooled domain6044 and586events',common['n']==rules['known_prior_common_n']==source['common_domain']['n']==6044 and
        common['all_cause_events8']==rules['known_prior_common_events8']==source['common_domain']['events8']==586 and
        eligible['n']==source['eligible_n'] and close(common['weight_sum'],source['common_domain']['weight_sum']))
    a.check('pooled sixyear weight convention not temporal32',r['weight_rule']=='(2/3)*WTSSCB4Y for1999/2001; (1/3)*WTSSCB2Y for2003')
    a.check('exact ten causes and separate unranked unknown',set(causes)==set(ALL_CODES) and unknown['clinically_ranked'] is False and unknown['support_evaluated'] is False)
    population_check(a,'common',common); population_check(a,'eligible',eligible)
    a.check('exact three baseline cycles',set(cycles)=={'1999','2001','2003'})
    for cycle in cycles:
        for domain in ('common','eligible'):population_check(a,cycle+' '+domain,cycles[cycle][domain])
    for domain,whole in [('common',common),('eligible',eligible)]:
        a.check(domain+' percycle population sums',all(sum(cycles[c][domain][k] for c in cycles)==whole[k] for k in
            ('n','all_cause_events8','alive_n','all_time_deaths_n','unknown_cause_deaths8','unknown_cause_deaths_after96')) and
            all(close(sum(cycles[c][domain][k] for c in cycles),whole[k]) for k in ('weighted_all_cause_events8','weight_sum')) and
            all(sum(cycles[c][domain]['all_time_cause_counts'][code] for c in cycles)==whole['all_time_cause_counts'][code] for code in ALL_CODES+('unknown',)))
    check_summary(a,'unknown overall',unknown,common,clinical=False)
    a.check('unknown early count equals common endpoint unknown',unknown['events8']==common['unknown_cause_deaths8'])
    for cycle in cycles:check_summary(a,'unknown '+cycle,unknown['per_cycle'][cycle],cycles[cycle]['common'],clinical=False)
    for code,group in causes.items():
        a.check(code+' exact id/ranking class',group['code']==code and group['cause_id']==definitions['groups'][code]['id'] and
            group['eligible_for_defined_clinical_ranking']==(code!='010'))
        check_summary(a,code+' overall',group['overall'],common,unknown_n=common['unknown_cause_deaths8'])
        for cycle in cycles:check_summary(a,code+' '+cycle,group['per_cycle'][cycle],cycles[cycle]['common'],unknown_n=cycles[cycle]['common']['unknown_cause_deaths8'])
        a.check(code+'8year cause counts no greater than alltime causes',group['overall']['events8']<=common['all_time_cause_counts'][code] and all(
            group['per_cycle'][cycle]['events8']<=cycles[cycle]['common']['all_time_cause_counts'][code] for cycle in cycles))
        a.check(code+' percycle event-weight moments add up',all(sum(group['per_cycle'][cycle][key] for cycle in cycles)==group['overall'][key] for key in
            ('n','events8','noncase_rows','label_unknown_rows','positive_weight_event_rows')) and
            all(close(sum(group['per_cycle'][cycle][key] for cycle in cycles),group['overall'][key]) for key in
                ('weighted_events8','event_weight_squared_sum','common_weight_sum')) and
            close(max((group['per_cycle'][cycle]['max_event_weight8'] for cycle in cycles if group['per_cycle'][cycle]['max_event_weight8'] is not None),default=None),group['overall']['max_event_weight8']))
    all_records=[g['overall'] for g in causes.values()]+[unknown]
    a.check('ten+unknown exact global event accounting',sum(x['events8'] for x in all_records)==common['all_cause_events8']==r['accounting']['ten_groups_plus_unknown_events8'] and
        close(sum(x['weighted_events8'] for x in all_records),common['weighted_all_cause_events8']) and
        close(r['accounting']['ten_groups_plus_unknown_weighted_events8'],common['weighted_all_cause_events8']) and r['accounting']['exact_raw_event_accounting'] is True)
    for cycle in cycles:
        records=[g['per_cycle'][cycle] for g in causes.values()]+[unknown['per_cycle'][cycle]]; population=cycles[cycle]['common']
        a.check(cycle+' ten+unknown count-weight accounting',sum(x['events8'] for x in records)==population['all_cause_events8'] and
            close(sum(x['weighted_events8'] for x in records),population['weighted_all_cause_events8']))
    # Independent outcome-free public masked-PSU mapping.
    inventory={str(s):[1,2,3] if s==1 else [1,2] for s in range(1,44)}
    a.check('full eligible common87PSU design preserved',r['design']==dict(strata_n=43,masked_PSU_n=87,eligible_common_full_inventory_equal=True,
        PSU_per_stratum_distribution={'2':42,'3':1}) and source['full_DEMO_masked_design']==source['eligible_masked_design']==inventory)
    keys=[s+':'+str(psu) for s,values in inventory.items() for psu in values]; outer=hashed_folds(keys,'b2m31:outer',5); folds=r['folds']
    a.check('independent outer mapping SHA and counts',folds['outer_mapping_sha256']==map_digest(outer) and
        folds['outer_evaluation_PSU_counts']=={str(f):sum(v==f for v in outer.values()) for f in range(5)} and
        folds['constructed_on_eligible_before_common'] is True and folds['outer_salt']=='b2m31:outer' and folds['inner_salt_pattern']=='b2m31:inner:{outer}')
    for record in folds['outer']:
        fold=record['outer_fold']; ikeys=[key for key,value in outer.items() if value!=fold]; inner=hashed_folds(ikeys,'b2m31:inner:'+str(fold),3)
        a.check('outer'+str(fold)+' independent inner mapping SHA and design',record['training_PSU_n']==len(inner) and
            record['evaluation_PSU_n']==sum(v==fold for v in outer.values()) and record['PSU_overlap_n']==0 and
            record['inner_mapping_sha256']==map_digest(inner) and record['inner_evaluation_PSU_counts']=={str(f):sum(v==f for v in inner.values()) for f in range(3)})
    expected_keys={(f,'outer',subset,None) for f in range(5) for subset in ('training','evaluation')}|{
        (f,'inner',subset,i) for f in range(5) for i in range(3) for subset in ('training','evaluation')}
    for code,group in causes.items():
        cells=group['fold_cells']; lookup={cell_key(c):c for c in cells}; overall=group['overall']
        a.check(code+' exact40fold cells',len(cells)==len(lookup)==40 and set(lookup)==expected_keys)
        a.check(code+' every cell distinguishes early unknown and known noncases',all(c['noncase_rows']==known_noncase(c['n'],c['events8'],c['label_unknown_rows']) and
            c['positive_weight_event_rows']==c['events8'] and 0<=c['event_strata_n']<=c['event_PSU_n']<=c['events8'] for c in cells))
        evals=[lookup[(f,'outer','evaluation',None)] for f in range(5)]
        additive=('n','events8','noncase_rows','label_unknown_rows','positive_weight_event_rows','event_PSU_n')
        wadditive=('weighted_events8','event_weight_squared_sum','common_weight_sum')
        a.check(code+' five outer evaluations reconstruct global support/moments',all(sum(c[k] for c in evals)==overall[k] for k in additive) and
            all(close(sum(c[k] for c in evals),overall[k]) for k in wadditive))
        for fold in range(5):
            ot=lookup[(fold,'outer','training',None)]; oe=lookup[(fold,'outer','evaluation',None)]
            ie=[lookup[(fold,'inner','evaluation',i)] for i in range(3)]
            a.check(code+' outer'+str(fold)+' disjoint train/eval and inner partitions',all(ot[k]+oe[k]==overall[k] for k in additive) and
                all(close(ot[k]+oe[k],overall[k]) for k in wadditive) and
                all(sum(c[k] for c in ie)==ot[k] for k in additive) and
                all(close(sum(c[k] for c in ie),ot[k]) for k in wadditive) and
                all(all(lookup[(fold,'inner','training',i)][k]+lookup[(fold,'inner','evaluation',i)][k]==ot[k] for k in additive) and
                    all(close(lookup[(fold,'inner','training',i)][k]+lookup[(fold,'inner','evaluation',i)][k],ot[k]) for k in wadditive) for i in range(3)))
        expected=independently_expected_failures(overall,cells)
        a.check(code+' all83 support comparisons and exact failed cells',group['support']['threshold_check_count']==83 and
            group['support']['thresholds_passed']==(not expected) and Counter(map(failure_key,group['support']['failures']))==Counter(map(failure_key,expected)))
    ranked=rank_defined({code:dict(weighted_events8=causes[code]['overall']['weighted_events8'],raw_events8=causes[code]['overall']['events8']) for code in DEFINED_CODES})
    a.check('independent9cause weighted ranking with exact tie rules',r['defined_cause_ranking']==ranked and all(causes[code]['defined_cause_rank']==i for i,code in enumerate(ranked,1)) and causes['010']['defined_cause_rank'] is None)
    a.check('fixed top5 selected before support and no010',r['priority_top5']==ranked[:5] and all(causes[code]['priority_top5']==(code in ranked[:5]) for code in ALL_CODES))
    unknown_block=common['unknown_cause_deaths8']>0
    expected_proceed=[code for code in ranked[:5] if causes[code]['support']['thresholds_passed'] and not unknown_block]
    a.check('top5 support filtering without backfill',r['proceed_codes']==expected_proceed and all(causes[code]['proceed_to_later_protocol']==(code in expected_proceed) for code in ALL_CODES))
    gate=r['gate']
    a.check('global known-cause gate separate from individual support',gate['source_structure_valid'] is True and gate['cause_labels_complete_for_common_horizon']==(not unknown_block) and
        gate['source_endpoint_qualified']==(not unknown_block) and gate['blockers']==(['unknown_common_cause_at_or_before96'] if unknown_block else []) and gate['no_comparative_fits_performed'] is True)
    a.check('later comparison availability matches at least one supportedpriority cause',gate['later_comparisons_allowed']==bool(expected_proceed))
    a.check('support blocker distinct from incomplete endpoint cause labels',gate['support_blockers']==(['no_supported_priority_cause'] if not unknown_block and not expected_proceed else []))
    # Completed scientific CI evidence supplied after completion; no network calls by this script.
    if Path(args.ci_jobs).exists():
        jobs=json.loads(Path(args.ci_jobs).read_text())
        a.check('final CI three jobs complete successfully',len(jobs)==3 and all(j['run_id']==args.run_id and j['status']=='completed' and j['conclusion']=='success' for j in jobs))
        for job in jobs:
            ownlog=log_text(Path(args.ci_log_dir)/(args.ci_log_prefix+str(job['id'])+'.log'))
            a.check(job['name']+' correct scientific checkout',args.science_commit in ownlog)
            if job['name'].startswith('synthetic'):
                a.check(job['name']+' testcount',bool(re.search(r'Ran '+str(args.tests_per_python)+r' tests in [\d.]+s\n\s*\nOK\n',ownlog)))
            else:a.check('two stdout blocks agree with independent completed empirical log',all(marked(ownlog,name)==marked(log,name) for name in paths))
    else:a.check('completed CI independent evidence supplied',False)
    def forbidden(value):
        if isinstance(value,dict):return [k for k in value if k.lower() in {'seqn','participant_ids','individual_labels','individual_predictions','participant_predictions','oof_predictions'}]+[k for v in value.values() for k in forbidden(v)]
        if isinstance(value,list):return [k for v in value for k in forbidden(v)]
        return []
    a.check('aggregate-only exports without participant or prediction keys',not forbidden(r) and not forbidden(manifest))
    a.summary=dict(science_commit=args.science_commit,plan_commit=PLAN_COMMIT,run_id=args.run_id,job_id=args.job_id,
        actual_PR_merge_checkout=(re.search(r'git log -1 --format=%H\n([0-9a-f]{40})',log).group(1)
            if re.search(r'git log -1 --format=%H\n([0-9a-f]{40})',log) else None),
        n=common['n'],allcause_events8=common['all_cause_events8'],unknown_early=common['unknown_cause_deaths8'],
        ranked_codes=ranked,priority_top5=ranked[:5],proceed_codes=expected_proceed,
        causes={code:dict(events8=causes[code]['overall']['events8'],event_PSU_n=causes[code]['overall']['event_PSU_n'],
            weighted_events8=causes[code]['overall']['weighted_events8'],weighted_incidence8=causes[code]['overall']['weighted_incidence8'],
            support_passed=causes[code]['support']['thresholds_passed'],failed_checks=len(causes[code]['support']['failures'])) for code in ALL_CODES},
        files_sha256={name:sha(path) for name,path in paths.items()},scientific_files=len(actual),historical_files=len(old),
        source_files=50,historical_data_files=25,tests_per_python=args.tests_per_python,job_log_sha256=sha(args.log))
    limits=[
        'This is a source/count/ranking/support gate only. No models, comparative prediction metrics, model bootstrap, Fine–Gray or joint integration are fitted or inferred.',
        'Source hashes, parsed outputs, accounting, ratios, fold maps, ranks and all support failures are independently reconciled from aggregates. No raw records or source byte downloads were read, and no second person-level cause reanalysis was performed.',
        'Top5 by weighted burden is fixed before support; failed priorities are not replaced. Count-based selection is explicit and does not justify confirmatory inference or exclusive marker-cause assignments.',
        'Unknown early causes block future comparisons and are never negative labels. Known competing deaths are noncases only for the fixed-horizon cumulative-incidence binary target, not censored or alive.',
        'Complete-case/LOD/linkage selection, reported-value assay/storage/calibration uncertainty and prior cohort use remain. Weights and pragmatic support thresholds do not certify power, EPV, clinical validity or national/global transport.',
        'Both historical mortality PDF denials remain denied_not_retried without new hashes or retries. Completed GitHub CI logs may be independently read; no raw/PDF/ZIP GET is used by this audit.'
    ]
    report=dict(schema='independent-cause33A-aggregate-audit/1',passed=sum(x['passed'] for x in a.checks),total=len(a.checks),
        all_passed=all(x['passed'] for x in a.checks),summary=a.summary,checks=a.checks,limits=limits)
    out=Path(args.out);out.with_suffix('.json').write_bytes(canonical(report))
    lines=['# Независимый аудит cause33A','',f"**{report['passed']}/{report['total']} PASS; общий исход: {'PASS' if report['all_passed'] else 'FAIL'}.**",'',
        f"Scientific commit `{args.science_commit}`; protocol `{PLAN_COMMIT}`; CIrun{args.run_id},job{args.job_id}.",
        f"Actual PR merge checkout `{a.summary['actual_PR_merge_checkout']}`; synthetic Python3.11/3.13: {args.tests_per_python} tests each, completed success.",'',
        f"Common n={common['n']},allcause8y={common['all_cause_events8']},unknown early={common['unknown_cause_deaths8']}.",
        'Weighted burden rank: '+', '.join(ranked)+'. Priority top5: '+', '.join(ranked[:5])+'. Proceed: '+(', '.join(expected_proceed) or 'none')+'.','',
        '| Code | Events8 | Event PSU | Weighted events8 | Support failures | Priority | Proceed |','|---|---:|---:|---:|---:|---|---|']
    for code in ALL_CODES:
        g=causes[code];o=g['overall']
        lines.append(f"| {code} | {o['events8']} | {o['event_PSU_n']} | {o['weighted_events8']:.8f} | {len(g['support']['failures'])} | {g['priority_top5']} | {g['proceed_to_later_protocol']} |")
    lines+=['','Проверены83supportcomparisons/cause,40foldcells/cause,all10+unknown accounting,percycleweights/ratios,independent hash folds,weightedranking/ties/no backfill,old95hashes иsource50/data25. Нет marker-benefit или clinical/external-validation вывода.','',
        'Ограничения:','']+['- '+s for s in limits]
    failed=[x['name'] for x in a.checks if not x['passed']]
    if failed:lines+=['','Failed checks:','']+['- '+s for s in failed]
    out.with_suffix('.md').write_text('\n'.join(lines)+'\n')
    print(json.dumps(dict(passed=report['passed'],total=report['total'],all_passed=report['all_passed'],failed=failed),ensure_ascii=False))
    return 0 if report['all_passed'] else 1

if __name__=='__main__':raise SystemExit(main())
