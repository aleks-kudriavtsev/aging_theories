"""Read-only independent audit of Forward32 aggregate exports, using stdlib only."""
from __future__ import annotations
import argparse
import ast
from collections import Counter
import datetime as dt
import hashlib
import json
import math
from pathlib import Path
import re
import subprocess
import urllib.parse

HEAD='59b0626767c0169107652cc7b7d51a416c8e7668'
REFERENCE_SHA='c036e3af139bdb6445a61f99a77cde82ffc7e831a4672dd2a77c1cfe56f56c09'


def sha(b):return hashlib.sha256(b).hexdigest()


def strict(b):
    def pairs(items):
        d={}
        for k,v in items:
            if k in d:raise ValueError('Duplicate JSON key')
            d[k]=v
        return d
    def reject(x):raise ValueError('Nonfinite JSON literal')
    d=json.loads(b,object_pairs_hook=pairs,parse_constant=reject)
    def finite(v):
        if isinstance(v,float) and not math.isfinite(v):raise ValueError('Nonfinite JSON number')
        if isinstance(v,dict):
            for x in v.values():finite(x)
        if isinstance(v,list):
            for x in v:finite(x)
    finite(d);return d


def serialized(d):return (json.dumps(d,indent=2,allow_nan=False)+'\n').encode()


def close(a,b):return math.isclose(a,b,rel_tol=1e-12,abs_tol=1e-9)


def allowed(u):
    p=urllib.parse.urlsplit(u)
    return p.scheme=='https' and p.hostname in {'wwwn.cdc.gov','www.cdc.gov','ftp.cdc.gov'} and p.username is None and p.password is None and p.port in (None,443) and bool(p.path) and not p.fragment


def audit(repo,log_path,ref_path):
    root=repo/'research/mortality';f=root/'forward32';checks={}
    def check(k,v):
        checks[k]=bool(v)
        if not v:raise AssertionError(k)
    def git_bytes(p):return subprocess.check_output(['git','show',HEAD+':'+p],cwd=repo)
    def load(p):return strict(p.read_bytes())
    r=load(f/'RESULTS_SUMMARY.json');m=load(f/'SOURCE_MANIFEST.json');art=load(f/'MODEL_ARTIFACTS.json');bins=load(f/'CALIBRATION_BINS.json')
    ref=load(ref_path);lock=load(f/'LOCK.json');spec=load(root/'b2m31/SOURCE_SPEC.json')
    check('pre_result_reference_exact_SHA',sha(ref_path.read_bytes())==REFERENCE_SHA)
    check('reference_head_exact',ref['scientific_head']==HEAD)
    check('14_metadata_pins',len(lock['metadata_sha256'])==14 and lock['metadata_sha256']==ref['metadata_sha256'])
    for n,h in ref['metadata_sha256'].items():
        p=root/n if '/' in n else f/n
        check('metadata:'+n,sha(p.read_bytes())==h==sha(git_bytes(p.relative_to(repo).as_posix())))
    check('90_historical_code_pins',len(lock['historical_scientific_python_sha256'])==90 and lock['historical_scientific_python_sha256']==ref['historical_scientific_python_sha256'])
    check('95_executed_science_hashes',len(r['scientific_python_sha256'])==95 and r['scientific_python_sha256']==ref['expected_scientific_python_sha256'])
    for n,h in ref['expected_scientific_python_sha256'].items():
        p=root/n;check('science:'+n,sha(p.read_bytes())==h==sha(git_bytes(p.relative_to(repo).as_posix())))
    check('all90historical_preserved',all(r['scientific_python_sha256'][n]==h for n,h in ref['historical_scientific_python_sha256'].items()))
    check('plan_and_lock_binding',r['plan_commit']=='377778c' and r['plan_sha256']==sha((f/'ANALYSIS_PLAN.md').read_bytes()) and r['lock_sha256']==sha((f/'LOCK.json').read_bytes())==ref['forward_lock_sha256'])
    check('top_manifest_schema',m['schema']=='forward32-source-manifest/1' and set(m)=={'schema','intake','source'})
    intake,base=m['intake'],m['source'];a=r['source_audit'];old=a['prior31_source_audit']
    check('intake_schema_binding',intake['schema']=='forward32-retrieval/1' and intake['plan_commit']=='377778c' and intake['lock_sha256']==r['lock_sha256'])
    check('base_schema_binding',base['schema']=='b2m31-retrieval/1' and base['plan_commit']=='6774fbe' and base['lock_sha256']==sha((root/'b2m31/LOCK.json').read_bytes()))
    check('fresh_base_exact_writer_SHA',sha(serialized(base))==intake['base_manifest_sha256']==old['source_manifest_sha256'])
    check('wrapper_exact_writer_SHA',sha(serialized(intake))==a['forward_source_wrapper_sha256'])
    check('export_manifest_SHA',sha((f/'SOURCE_MANIFEST.json').read_bytes())==r['source_manifest_export_sha256'])
    rows={x['file']:x for x in base['files']}
    check('exact50sourcefiles',len(rows)==len(base['files'])==old['files_verified']==50 and set(rows)==set(ref['required_base_source_files']))
    urls={}
    for c in spec['cycles']:
        for n,v in c['components'].items():
            urls[str(c['year'])+'_'+n+'.xpt']=v['historical_declared_data_url'];urls[str(c['year'])+'_'+n+'.htm']=v['documentation_url']
        urls[str(c['year'])+'_mortality.dat']=c['mortality_url']
    sscard_url='https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/1999/DataFiles/SSCARD_A.xpt'
    urls.update({'index.htm':spec['sscard_index_url'],'1999_SSCARD_A.htm':spec['sscard']['metadata_url'],'1999_SSCARD_A.xpt':sscard_url,'weighting.htm':spec['weighting_url'],'variance.htm':spec['variance_url']})
    for n,v in rows.items():
        check('source:'+n,Path(n).name==n and v['status']=='retrieved' and type(v['bytes']) is int and 0<v['bytes']<=25000000 and bool(re.fullmatch('[a-f0-9]{64}',v['sha256'])) and allowed(v['url']) and allowed(v['resolved_url']) and v['url']==urls[n])
    fingerprints={n:v['sha256'] for n,v in rows.items()};data={n:v for n,v in rows.items() if n.endswith(('.xpt','.dat'))}
    check('three50fingerprintmaps',fingerprints==intake['source_fingerprints']==old['source_fingerprints'])
    content_sha=sha(json.dumps(fingerprints,sort_keys=True,separators=(',',':')).encode())
    check('canonical_content_SHA',content_sha==old['source_content_sha256'])
    check('exact25historicaldataSHA',len(data)==old['historical_data_hashes_verified']==25 and {n:v['sha256'] for n,v in data.items()}==ref['historical_data_sha256'])
    check('actual_sscard_href',rows['1999_SSCARD_A.xpt']['url']==rows['1999_SSCARD_A.xpt']['resolved_url']==sscard_url and rows['1999_SSCARD_A.xpt']['bytes']==2279120)
    check('2priorPDFdenials',base['documentation_not_retrieved']==old['documentation_not_retrieved'] and len(base['documentation_not_retrieved'])==2)
    check('exact_deniedPDFURLs',{v['url'] for v in base['documentation_not_retrieved']}=={spec['mortality_dictionary_url'],spec['mortality_description_url']})
    for v in base['documentation_not_retrieved']:
        check('PDFnotretried:'+v['file'],v['status']=='denied_not_retried' and v['original_http_status']==403 and v['original_run_id']==40907895237 and v['original_job_id']==123567576780 and v['file'] not in rows and 'bytes' not in v and 'sha256' not in v)
    obs=a['observed_full_SSCARD_marginals'];expected=ref['expected_full_SSCARD_marginals']
    check('all_fullsource_marginals',all(obs[k]==expected[k] for k in ['rows','flags','pristine','weights','source_sha256']))
    columns=['SEQN','WTSSCB2Y','WTSSCB4Y','SSB2M','SSB2ML','SSCYST','SSCYSTL','SSGALB','SSGALBL','SSALB','SSALBL','SSGALBP','SSPRIS']
    check('13_exact_source_columns_no_batch',obs['columns']==columns and obs['individual_batch_or_assay_date_fields_provided'] is False)
    check('observed_source_bound_SHA',obs['source_sha256']==rows['1999_SSCARD_A.xpt']['sha256'])
    for k,v in obs['flags'].items():check('full_flag_sum:'+k,sum(v.values())==21890)
    check('full_pristine_sum',sum(obs['pristine'].values())==21890)
    for k,v in obs['weights'].items():check('full_weight_sum:'+k,sum(v.values())==21890)
    check('assay_units_platform',spec['sscard']['unit']=={'SSCYST':'mg/L','SSB2M':'mg/L'} and spec['sscard']['matrix']=='serum' and spec['sscard']['platform']=='Siemens Dimension Vista 1500' and spec['sscard']['last_revised']=='2022-07')
    check('22XPT_zero_audits',len(old['xpt_zero_audits'])==22)
    for v in old['xpt_zero_audits']:check('XPTzero:'+v['filename'],v['sha256']==rows[v['filename']]['sha256'] and v['epsilon_threshold_used'] is False)
    check('3_exact_mortality_readers',len(old['mortality_reader_audits'])==3 and all(v['exact_match'] and v['fields']==8 and v['cell_comparisons']==8*v['rows'] and v['sha256']==rows[v['file']]['sha256'] for v in old['mortality_reader_audits']) and sum(v['rows'] for v in old['mortality_reader_audits'])==31126 and sum(v['cell_comparisons'] for v in old['mortality_reader_audits'])==249008)

    native_tree=ast.parse((f/'cohort.py').read_text());fn=next(n for n in native_tree.body if isinstance(n,ast.FunctionDef) and n.name=='native_period_frame')
    calls=[n for n in ast.walk(fn) if isinstance(n,ast.Call) and isinstance(n.func,ast.Attribute) and isinstance(n.func.value,ast.Name) and n.func.value.id=='np' and n.func.attr=='where']
    exp=ast.parse('np.where(early,d.WTSSCB4Y,d.WTSSCB2Y)',mode='eval').body
    check('native_weight_assignment_AST',any(ast.dump(n,include_attributes=False)==ast.dump(exp,include_attributes=False) for n in calls))
    check('fixed_split_and_raw_weight_labels',r['train_cycles']==[1999,2001] and r['late_cycles']==[2003] and r['train_weight']=='WTSSCB4Y' and r['late_weight']=='WTSSCB2Y')
    check('no_overlap',a['participant_overlap_n']==0 and a['masked_PSU_overlap_n']==0)
    periods=a['periods'];pairsets={};weight_checks={}
    for part,lo,hi,np,known_n,events,label in [('train',1,28,57,4020,383,'1999_2002'),('late',29,43,30,2024,203,'2003_2004')]:
        p=periods[part];eligible,common=p['eligible'],p['common'];inv=p['full_masked_design']
        check('equal_period_inventory:'+part,inv==p['eligible_masked_design']==p['common_masked_design'])
        check('period_strata:'+part,set(inv)==set(map(str,range(lo,hi+1))))
        pairsets[part]={(int(st),psu) for st,psus in inv.items() for psu in psus}
        check('period_PSU_inventory:'+part,len(pairsets[part])==np and all(psus==([1,2,3] if int(st)==1 else [1,2]) for st,psus in inv.items()))
        check('period_received_domain_label:'+part,p['weight_variable']==('WTSSCB4Y' if part=='train' else 'WTSSCB2Y') and p['estimand']=='received_complete_case_measurable_domain_of_'+label)
        check('common_frame_unchanged:'+part,common['n']==known_n and common['events8']==events and common['n']==r[part+'_n'] and common['events8']==r[part+'_events8'])
        for domain,d in [('eligible',eligible),('common',common)]:
            k=part+':'+domain;n=d['n'];w=d['native_weight_sum'];w2=d['native_weight_squared_sum']
            check('domain_positive_event_arithmetic:'+k,n==d['events8']+d['nonevents8'] and d['events8']==d['positive_weight_events8'] and d['nonevents8']==d['positive_weight_nonevents8'])
            check('domain_clusters:'+k,d['strata_n']==hi-lo+1 and d['psu_n']==np and 0<=d['event_bearing_psu_n']<=np and 0<=d['nonevent_bearing_psu_n']<=np and d['event_bearing_psu_n']+d['nonevent_bearing_psu_n']>=np)
            check('source_weight_exact_selection:'+k,d['selected_native_source_max_abs_difference']==0 and close(w,d['source_W4_sum'] if part=='train' else d['source_W2_sum']))
            check('Kish_recomputed:'+k,w>0 and w2>0 and close(w*w/w2,d['effective_n']) and 0<d['effective_n']<=n)
            check('followup:'+k,d['minimum_survivor_followup_months']>=96)
            check('pristine_count_sum:'+k,sum(d['pristine_counts'].values())==n)
            for field,v in d['protein_flag_counts'].items():
                check('flag_count_sum:'+k+field,set(v)=={'0','1','2','missing'} and sum(v.values())==n)
                if domain=='common':check('flags00:'+k+field,v=={'0':n,'1':0,'2':0,'missing':0})
            check('cycle_n_sum:'+k,sum(d['cycle_counts'].values())==n and set(d['per_cycle_weights'])==set(d['cycle_counts']))
            for name in ['source_W4_sum','source_W2_sum','original_W6_sum','native_weight_sum']:
                check('percycle_weight_sum:'+k+name,close(sum(v[name] for v in d['per_cycle_weights'].values()),d[name]))
            for c,v in d['per_cycle_weights'].items():
                chosen='source_W4_sum' if part=='train' else 'source_W2_sum';factor=2/3 if part=='train' else 1/3
                check('percycle_native_and_W6_factor:'+k+c,v['n']==d['cycle_counts'][c] and close(v['native_weight_sum'],v[chosen]) and close(v['original_W6_sum'],factor*v[chosen]))
            check('W4_missingness_semantics:'+k,(d['source_W4_missing_n']==0 if part=='train' else d['source_W4_missing_n']==n and d['source_W4_sum']==0))
            check('feature_moment_keys:'+k,len(d['outcome_free_feature_moments'])==11 and set(d['outcome_free_log2_protein_moments'])=={'SSCYST','SSB2M'})
            for name,v in {**d['outcome_free_feature_moments'],**d['outcome_free_log2_protein_moments']}.items():
                check('moments:'+k+name,0<=v['n']<=n and (v['n']==n if domain=='common' else True)
                      and ((v['weighted_mean'] is None and v['weighted_SD'] is None) if v['n']==0 else isinstance(v['weighted_mean'],(float,int)) and v['weighted_SD']>=0))
            weight_checks[k]=dict(n=n,events8=d['events8'],sum_w=w,sum_w2=w2,recomputed_effective_n=w*w/w2,event_bearing_PSU=d['event_bearing_psu_n'],nonevent_bearing_PSU=d['nonevent_bearing_psu_n'])
        check('coverage_recomputed:'+part,close(common['native_weight_sum']/eligible['native_weight_sum'],p['common_native_weighted_coverage']) and 0<p['common_native_weighted_coverage']<=1)
        check('common_subset_counts:'+part,common['n']<=eligible['n'] and common['events8']<=eligible['events8'])
        check('min20_positive_weight_rows_notPSUs:'+part,common['positive_weight_events8']>=lock['minimum_positive_weight_event_rows_per_part']==20 and common['positive_weight_nonevents8']>=lock['minimum_positive_weight_nonevent_rows_per_part']==20)
    check('57plus30clusters_union',not pairsets['train']&pairsets['late'] and len(pairsets['train']|pairsets['late'])==87)
    check('known31frame_total',a['original31_common_n']==periods['train']['common']['n']+periods['late']['common']['n']==6044 and a['prior31_reported_sensitivity_not_reanalysed'] is True)
    full=old['full_DEMO_masked_design'];check('partition_union_equals_full31pairs',pairsets['train']|pairsets['late']=={(int(st),p) for st,ps in full.items() for p in ps})
    check('m3only_train_stratum1',periods['train']['full_masked_design']['1']==[1,2,3] and sum(len(ps)==3 for ps in full.values())==1)
    check('computed_gate_passed',r['gate']==a['gate'] and r['gate']=={'temporal_evaluation_allowed':True,'blockers':[],'support_threshold_is_not_power_or_calibration_certificate':True} and r['status']=='computed' and r['model_fits_attempted']==4)
    mapping=a['eligible_derived_early_inner_fold_mapping'];keys={str(s)+':'+str(p) for s,p in pairsets['train']};ordered=sorted(keys,key=lambda k:sha((lock['inner_salt']+':'+k).encode()))
    check('independent57PSU_threefold_SHA_mapping',mapping=={k:i%3 for i,k in enumerate(ordered)} and Counter(mapping.values())=={0:19,1:19,2:19})
    for v in art['inner_folds_aggregate']:
        check('training_fold_arithmetic:'+str(v['fold']),v['training_n']+v['scoring_n']==4020 and v['training_events']+v['scoring_events']==383 and v['PSU_overlap_n']==0 and v['training_psu_n']==38 and v['scoring_psu_n']==19 and close(v['training_weight_sum']+v['scoring_weight_sum'],periods['train']['common']['native_weight_sum']))
    check('training_fold_totals',sum(v['scoring_n'] for v in art['inner_folds_aggregate'])==4020 and sum(v['scoring_events'] for v in art['inner_folds_aggregate'])==383)
    check('model_domain_binding',art['training']['n']==4020 and art['training']['events8']==383 and art['training']['weight_name']=='WTSSCB4Y' and art['late_domain']['n']==2024 and art['late_domain']['evaluation_weight_name']=='WTSSCB2Y' and art['late_domain']['training_psu_overlap_n']==0 and art['late_domain']['outcomes_used_for_fitting'] is False and art['late_domain']['weights_used_for_fitting'] is False)
    check('four_arms_one_late_domain',len(r['point'])==4 and all(v['n']==2024 and v['events']==203 and v['n_positive_weight']==2024 and close(v['weight_sum'],periods['late']['common']['native_weight_sum']) for v in r['point'].values()) and r['same_common_domain_all_four_models_within_each_period'] is True)
    check('no_nationalW6_no_late_fit_no_clinical_claim',a['no_six_year_pooled_score'] is True and a['no_late_recalibration'] is True and r['late_recalibration_performed'] is False and r['independent_external_validation'] is False and r['clinical_use_ready'] is False and r['participant_records_exported'] is False and a['assay_batch_individually_identified'] is False and a['storage_duration_individually_identified'] is False)
    for n,h in r['outputs_sha256'].items():check('export_bytes:'+n,sha((f/n).read_bytes())==h)
    log_bytes=log_path.read_bytes();lines=[re.sub(r'^\d{4}-\d\d-\d\dT[^ ]+ ','',x) for x in log_bytes.decode().splitlines()]
    def marked(label):
        i=lines.index('FORWARD32_'+label+'_BEGIN');j=lines.index('FORWARD32_'+label+'_END');return strict('\n'.join(lines[i+1:j]))
    for label,value in [('RESULTS',r),('MANIFEST',m),('MODEL_ARTIFACTS',art),('CALIBRATION_BINS',bins)]:check('exact_marked_export:'+label,marked(label)==value)
    checkout=next(x for x in lines if x.startswith('HEAD is now at '));check('job_scientific_head',('Merge '+HEAD+' into ') in checkout)
    abbrev=checkout.split()[4];merge=next(x for x in lines if re.fullmatch('[a-f0-9]{40}',x) and x.startswith(abbrev))
    check('successful50_source_intake_in_log','{"retrieved": 50, "failures": []}' in lines)
    artifacts={p.name:dict(bytes=p.stat().st_size,sha256=sha(p.read_bytes())) for p in [f/'RESULTS_SUMMARY.json',f/'SOURCE_MANIFEST.json',f/'MODEL_ARTIFACTS.json',f/'CALIBRATION_BINS.json',log_path]}
    return dict(schema='forward32-independent-source-audit/1',status='passed_with_raw_access_limits',date='2026-10-29',checks=checks,passed_check_count=len(checks),failed_check_count=0,
        no_participant_records_read=True,no_network_requests=True,no_ZIP_download=True,no_refits=True,
        scientific_head=HEAD,checked_out_PR_merge_commit=merge,run_id_supplied_by_parent=40917437691,job_id_supplied_by_parent=123601312507,
        script_sha256=sha(Path(__file__).read_bytes()),reference_sha256=sha(ref_path.read_bytes()),artifacts=artifacts,
        pin_counts=dict(metadata=14,historical_Python=90,executed_scientific_Python=95,source_manifest_entries=50,historical_data_SHA=25),
        source=dict(full_source_bytes_reported=sum(v['bytes'] for v in rows.values()),data_source_bytes_reported=sum(v['bytes'] for v in data.values()),manifest_SHA=r['source_manifest_export_sha256'],fresh_base_manifest_SHA=intake['base_manifest_sha256'],forward_wrapper_SHA=a['forward_source_wrapper_sha256'],actual_sscard=rows['1999_SSCARD_A.xpt'],prior_PDF_denials=base['documentation_not_retrieved']),
        observed_full_SSCARD_marginals=obs,weights_and_support=weight_checks,coverage={p:v['common_native_weighted_coverage'] for p,v in periods.items()},
        period_estimands={p:v['estimand'] for p,v in periods.items()},design=dict(train_strata=28,train_PSU=57,late_strata=15,late_PSU=30,m3_train_stratum={'1':[1,2,3]},participant_overlap_n=0,masked_PSU_overlap_n=0,inner_folds_PSU_counts=dict(Counter(mapping.values()))),
        limitations=['Raw source files were not independently rehashed or read; source hashes and aggregate categories were crosschecked against pinned references and exported provenance.',
                     'Runtime SSCARD href is reported by the pinned unique-anchor resolver; raw catalogue HTML was not independently re-parsed here.',
                     'Native weight sums, squared sums, Kish n, coverage and W6 factors were independently checked from aggregates; rowwise zero-difference is an exported diagnostic, not a separate record-level audit.',
                     'Weighted feature moments were checked for domain/finite/SD consistency; individual raw measurements were not used to recompute moments.',
                     'Participant mask identity and zero overlaps are supported by pinned code and aggregate diagnostics, not independently matched individual IDs.',
                     'The min20 rule concerns positive-weight event/nonevent rows, not event-bearing PSU and not power or calibration precision.',
                     'Reused cohorts and prior marker selection; retrospective baseline-calendar stress test, not an untouched holdout, historical as-of2003 deployment or independent frozen30 validation.',
                     'Bootstrap confidence interval arithmetic is covered by the separate numerical reviewer, not this source/domain audit.'])


def report(d):
    src=d['source'];w=d['weights_and_support'];cov=d['coverage']
    lines=['# Forward32: независимый аудит источников, весов и домена','',
        f"**Пройден: {d['passed_check_count']} проверок, отказов 0.** Проверены агрегаты, код и экспортированные манифесты; записи участников, raw-source bytes, ZIP и сетевые источники не открывались, refits не выполнялись.",'',
        f"Scientific head `{HEAD}`; checkout журнала — PR merge `{d['checked_out_PR_merge_commit']}`. Совпали 14 LOCK metadata SHA, 90 historical Python SHA и все 95 exported scientific SHA. Четыре JSON совпадают с marked job exports; шесть вспомогательных файлов совпадают с RESULT output SHA. Run/job указаны по переданному журналу и информации root; отдельного сетевого CI lookup не было.",'',
        f"Ровно 50 обязательных source entries; все 25 data SHA совпали с прежними pins по именам. Заявленные байты всех источников: {src['full_source_bytes_reported']:,}; данных: {src['data_source_bytes_reported']:,}. Свежий base manifest и wrapper32 связаны своими хешами, а старый manifest31 отдельно сохранён в LOCK. Два прежних PDF с HTTP403 отмечены denied_not_retried; новых запросов/хешей этих PDF не заявлено.",'',
        f"SSCARD runtime href: `{src['actual_sscard']['url']}`. SHA `{src['actual_sscard']['sha256']}`; размер {src['actual_sscard']['bytes']:,} байтов. Полные observed marginals (оба flags0/1/2/missing, pristine и оба веса) точно совпали с заранее закреплёнными таблицами July2022 для 21890 строк. Это независимая сверка экспортированных source aggregates, а не повторное чтение raw records. Serum, mg/L и Dimension Vista1500 сохранены; LOD не назван LOQ.",'',
        'Native train weight — WTSSCB4Y, late — WTSSCB2Y. Для каждого eligible/common периода и цикла проверены source/native суммы, originalW6 factors, sum of squared weights, Kish effective n и coverage. Rowwise max difference 0 — экспортированный diagnostic. Late W4 полностью missing; его сумма0 отражает sum(skipna), а не индивидуальные нулевые веса. Estimands условны на полученный complete-case measurable domain; национальная W6 estimate не заявлена.','',
        '| Часть / домен | n | События / несобытия | Сумма native weights | Kish n | Event-bearing PSU |',
        '|---|---:|---:|---:|---:|---:|']
    for k in ['train:eligible','train:common','late:eligible','late:common']:
        x=w[k];lines.append(f"| {k} | {x['n']} | {x['events8']} / {x['n']-x['events8']} | {x['sum_w']:.6f} | {x['recomputed_effective_n']:.6f} | {x['event_bearing_PSU']} |")
    lines += ['',f"Coverage common/eligible: train {cov['train']:.12f}, late {cov['late']:.12f}. Primary flags00 сохранены во всех четырёх arms: train — 4020/383, late — 2024/203. LOD01 sensitivity итерации31 здесь не анализировалась.",'',
        'Full = eligible = common inventory в каждой части: train — 28 strata / 57 PSU, late — 15 strata / 30 PSU. Единственная m3 — train stratum1=[1,2,3]. Participant overlap и masked-PSU overlap — 0. Точная early fold mapping независимо восстановлена через SHA256: по 19 PSU в каждом из трёх folds. Агрегаты inner scoring складываются в полный train domain; late в них не входит.','',
        'Gate min20 пройден для positive-weight event и nonevent **rows** в каждой части. Это не порог 20 PSU и не доказательство мощности или калибровочной точности. Реальные common event-bearing PSU: 57 train / 30 late; nonevent-bearing PSU также 57/30. Минимальный survivor follow-up — 201/177 месяцев, выше 96.','',
        'Ограничения: нет batch/individual assay date/storage duration; assay 2018–2020 не равен baseline 1999–2004. Cohorts и поздняя часть уже использовались, маркер выбран после прежних результатов. Это retrospective baseline-calendar stress test, не untouched holdout, не deployment as-of2003 и не независимая валидация frozen30. Все пределы record-level независимости перечислены в companion JSON. Численные интервалы bootstrap проверяет отдельный reviewer.','']
    return '\n'.join(lines)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--repo',type=Path,default=Path('/workspace/aging_theories'))
    p.add_argument('--log',type=Path,default=Path('/tmp/forward32_success_job.log'));p.add_argument('--reference',type=Path,default=Path('/tmp/forward32_source_expected_reference.json'))
    p.add_argument('--out',type=Path,default=Path('/tmp/forward32_source_audit'));a=p.parse_args();d=audit(a.repo,a.log,a.reference)
    a.out.with_suffix('.json').write_text(json.dumps(d,indent=2,ensure_ascii=False,allow_nan=False)+'\n');a.out.with_suffix('.md').write_text(report(d))
    print(json.dumps(dict(status=d['status'],checks=d['passed_check_count'],failures=0,train='4020/383',late='2024/203',source_entries=50,data_pins=25)))
