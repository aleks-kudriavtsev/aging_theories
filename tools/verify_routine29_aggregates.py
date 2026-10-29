"""Independent stdlib-only audit of routine29 exported aggregates; never fits models."""
from pathlib import Path
from collections import Counter
import argparse, csv, hashlib, json, math, re, subprocess

_parser=argparse.ArgumentParser(description=__doc__)
_parser.add_argument('--repo',type=Path,default=Path(__file__).resolve().parents[1])
_parser.add_argument('--job-log',type=Path,required=True)
_parser.add_argument('--out',type=Path,required=True)
_args=_parser.parse_args()
REPO=_args.repo.resolve()
OUT=_args.out
if OUT.exists():raise ValueError('New audit output directory required')
OUT.mkdir(parents=True)
PARENT=REPO/'research/mortality'
ROOT=PARENT/'routine29'
LOG=_args.job_log
RUN=40897873302
JOB=123531809031
ARMS={'simple','simple_all3','routine','routine_all3','rank_all3','clip_all3'}
PRIMARY={'routine','routine_all3'}
checks={}; details={}; errors=[]

def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def check(name,condition,detail=None):
    checks[name]=bool(condition)
    if detail is not None:details[name]=detail
    if not condition:errors.append(name)

def same(a,b,tolerance=1e-12):
    return a is not None and b is not None and math.isfinite(float(a)) and math.isfinite(float(b)) and abs(float(a)-float(b))<=tolerance

def marker(text,label):
    lines=text.splitlines(); begin='ROUTINE29_'+label+'_BEGIN';end='ROUTINE29_'+label+'_END'
    if lines.count(begin)!=1 or lines.count(end)!=1:raise ValueError('Nonunique block: '+label)
    a,b=lines.index(begin),lines.index(end)
    if a>=b:raise ValueError('Block order: '+label)
    return '\n'.join(lines[a+1:b])+'\n'

def percentile(values,p):
    # Independent order-statistic interpolation, without numpy or pipeline imports.
    ordered=sorted(values);position=(len(ordered)-1)*p
    lower=math.floor(position);upper=math.ceil(position)
    return ordered[lower]+(ordered[upper]-ordered[lower])*(position-lower)

def metric_checks(prefix,m,n,events):
    check(prefix+'_n_events',m['n']==n and m['events']==events)
    check(prefix+'_bias_arithmetic',same(m['bias_pp'],100*(m['predicted']-m['observed'])))
    check(prefix+'_OE_arithmetic',same(m['observed_expected'],m['observed']/m['predicted']))
    check(prefix+'_ranges',all(math.isfinite(m[k]) and 0<=m[k]<=1 for k in ('auc','brier','observed','predicted')))
    status=m['calibration_status']
    check(prefix+'_calibration_status',status!='computed' or (all(m[k] is not None and math.isfinite(m[k]) for k in ('calibration_in_large','calibration_intercept','calibration_slope')) and m['calibration_gradient_max']<=1e-5))

def run():
    plain='\n'.join(re.sub(r'^\S+ ','',line) for line in LOG.read_text().splitlines())
    files={'RESULTS':'RESULTS_SUMMARY.json','MANIFEST':'SOURCE_MANIFEST.json','PARAMETERS':'FOLD_PARAMETERS.json','CALIBRATION_BINS':'CALIBRATION_BINS.json'}
    decoded={}
    for label,name in files.items():
        value=json.loads(marker(plain,label));local=json.loads((ROOT/name).read_text())
        check(label.lower()+'_matches_independently_fetched_successful_job',local==value)
        canonical=json.dumps(value,ensure_ascii=False,indent=2,allow_nan=False)+'\n'
        check(label.lower()+'_bytes_equal_canonical_stdout',(ROOT/name).read_bytes()==canonical.encode())
        decoded[label]=local
    result=decoded['RESULTS'];manifest=decoded['MANIFEST'];fits=decoded['PARAMETERS'];bins=decoded['CALIBRATION_BINS']
    csv_stdout=marker(plain,'PRIMARY_BOOTSTRAP_CSV');local_csv=(ROOT/'primary_bootstrap.csv').read_text()
    check('primary_csv_bytes_equal_successful_stdout',local_csv==csv_stdout)
    check('primary_csv_hash_matches_result',hashlib.sha256(local_csv.encode()).hexdigest()==result['primary_bootstrap_csv_sha256'])
    records=list(csv.DictReader(local_csv.splitlines()))
    check('500_unique_prespecified_replicates',len(records)==500 and sorted(int(r['replicate']) for r in records)==list(range(500)))
    check('replicate_statuses_explicit',all(r['status'] in ('computed','undefined') for r in records))
    failures=sum(r['status']=='undefined' for r in records)
    check('replicate_failures_match_result',failures==result['bootstrap_failures'])
    independent_limits={}
    for key in ('delta_brier','delta_auc'):
        values=[float(r[key]) for r in records if r[key].strip() and math.isfinite(float(r[key]))]
        saved=result['primary_intervals'][key]
        check(key+'_defined_undefined_counts',len(values)==saved['valid'] and 500-len(values)==saved['undefined'])
        limits=[percentile(values,.025),percentile(values,.975)] if values else [None,None]
        independent_limits[key]=limits
        check(key+'_independent_linear_percentiles',all(same(x,y,1e-14) for x,y in zip(limits,[saved['lower'],saved['upper']])) if values else saved['lower'] is None and saved['upper'] is None,limits)
        check(key+'_defined_only_label',saved['defined_draws_only']==(len(values)!=500))
        check(key+'_fixed_partition_refit_scope',saved['scope']=='nested_preprocessing_tuning_refitting_with_fixed_partition')
    check('replicate_status_and_missing_contrasts_agree',all((r['status']=='computed' and r['delta_brier'].strip() and r['delta_auc'].strip()) or (r['status']=='undefined' and not r['delta_brier'].strip() and not r['delta_auc'].strip()) for r in records))
    check('all_interval_counts_labels_consistent',all(v['valid']+v['undefined']==500 and v['defined_draws_only']==(v['undefined']>0) for v in result['primary_intervals'].values()))
    check('bootstrap_seed_count_locked',result['bootstrap_seed']==20261101 and result['bootstrap_replicates']==500)

    scientific=result['scientific_python_sha256']
    check('all_recorded_scientific_hashes_match_current_tree',all(digest(PARENT/name)==value for name,value in scientific.items()))
    check('all_scientific_python_files_recorded',set(scientific)=={str(p.relative_to(PARENT)) for p in PARENT.rglob('*.py')})
    lock=json.loads((ROOT/'LOCK.json').read_text())
    check('lock_hash_matches_result',digest(ROOT/'LOCK.json')==result['lock_sha256'])
    check('plan_hash_matches_result',digest(ROOT/'ANALYSIS_PLAN.md')==result['plan_sha256'])
    check('metadata_hashes_match_lock',all(digest(ROOT/name)==value for name,value in lock['study_metadata_sha256'].items()))
    check('prior28_manifest_hash_locked',digest(PARENT/'binary28/SOURCE_MANIFEST.json')==lock['prior28_source_manifest_sha256'])
    frozen=json.loads((PARENT/'temporal26/LOCK.json').read_text())['models']
    check('historical_frozen_models_preserved',all(digest(PARENT/name)==value for name,value in frozen.items()))
    oldnames=subprocess.check_output(['git','ls-tree','-r','--name-only',lock['base_commit'],'--','research/mortality/binary28'],cwd=REPO,text=True).splitlines()
    check('all_binary28_files_preserved_from_locked_base',all((REPO/name).read_bytes()==subprocess.check_output(['git','show',lock['base_commit']+':'+name],cwd=REPO) for name in oldnames),{'files':len(oldnames)})
    check('plan_commit_ancestor_of_scientific_head',subprocess.run(['git','merge-base','--is-ancestor','7800458','b8d05a0'],cwd=REPO,capture_output=True).returncode==0)

    jobs=json.loads(Path('/tmp/routine29_jobs.json').read_text())['jobs']
    check('all_three_scientific_CI_jobs_success',len(jobs)==3 and all(j['status']=='completed' and j['conclusion']=='success' for j in jobs))
    for version,file in [('3.11','/tmp/routine29_synthetic311.log'),('3.13','/tmp/routine29_synthetic313.log')]:
        log=Path(file).read_text()
        check('CI_'+version+'_53_tests_pass',re.search(r'Ran 53 tests in',log) is not None and re.search(r'(?m)^\S+ OK$',log) is not None)

    audit=result['source_audit'];fingerprints={r['file']:r['sha256'] for r in manifest['files']}
    check('routine_manifest_hash_matches_result',digest(ROOT/'SOURCE_MANIFEST.json')==audit['routine_source_manifest_sha256'])
    check('manifest_unique_sources',len(fingerprints)==len(manifest['files']))
    check('all_sources_retrieved_public_CDC',all(r['status']=='retrieved' and r['url'].startswith(('https://wwwn.cdc.gov/','https://ftp.cdc.gov/')) and r['resolved_url'].startswith(('https://wwwn.cdc.gov/','https://ftp.cdc.gov/')) for r in manifest['files']))
    check('routine_source_fingerprints_equal_manifest',fingerprints==audit['routine_source_fingerprints'])
    stable=hashlib.sha256(json.dumps(fingerprints,sort_keys=True,separators=(',',':')).encode()).hexdigest()
    check('routine_stable_source_digest_matches',stable==audit['routine_source_content_sha256'])
    check('old_source_fingerprints_preserved',all(fingerprints[name]==value for name,value in audit['source_fingerprints'].items()))
    prior=json.loads((PARENT/'binary28/SOURCE_MANIFEST.json').read_text())
    check('prior28_scientific_source_bytes_preserved',all(fingerprints[r['file']]==r['sha256'] for r in prior['files'] if r['file']!='index.htm'))
    for kind,n,fields,cells,end,width in [('lab',29314,6,175884,1970,1977),('exam',31311,4,125244,1533,6233)]:
        a=audit['routine_readers'][kind]
        check(kind+'_independent_reader_exact_cells',a['exact_match'] is True and a['rows']==n and a['fields']==fields and a['cell_comparisons']==cells and a['expected_rows']==n)
        check(kind+'_source_and_SAS_hashes_match_manifest',a['source_sha256']==fingerprints[kind+'.dat'] and a['sas_sha256']==fingerprints[kind+'.sas'])
        check(kind+'_width_contract_matches_locked_metadata',a['SAS_MISSOVER'] is True and a['selected_fields_end']==end and a['nominal_content_width']==width and a['FORMAT_scaling_applied'] is False)
        widths=a['source_record_width_counts']
        check(kind+'_width_counts_preserve_all_rows',sum(widths.values())==n and all(end<=int(k)<=width for k in widths))
        # LRECL in the SAS documentation includes CRLF, but actual released
        # byte counts need not. Here the reported bytes reconcile to one byte
        # per physical record beyond the measured content widths.
        content_bytes=sum(int(k)*v for k,v in widths.items())
        check(kind+'_source_bytes_reconcile_one_byte_line_end',content_bytes+n==a['source_bytes'],{'content_bytes':content_bytes,'source_bytes':a['source_bytes'],'bytes_per_record_beyond_content':(a['source_bytes']-content_bytes)/n})
        check(kind+'_short_unselected_suffix_count_reconciles',sum(v for k,v in widths.items() if int(k)<width)==a['rows_with_absent_trailing_unselected_characters'])
    check('clinical_reader_exact_200500_cells',audit['clinical_reader']['rows']==20050 and audit['clinical_reader']['fields']==10 and audit['clinical_reader']['cell_comparisons']==200500 and audit['clinical_reader']['exact_match'])
    check('qualification_flow_preserved',audit['eligible_protein_weight_n']==3899 and audit['qualified_protein_n']==3864 and audit['common_complete_n']==3857 and audit['old28_common_n']==3857 and audit['old28_common_events8']==612)
    patterns=audit['exclusion_reason_patterns'];oldpatterns=audit['exclusion_patterns_in_old28_common']
    check('main_flow_3899_to_common_reconciles',sum(patterns.values())==3899 and patterns['complete']==result['n']==audit['routine_common_n'])
    check('old28_subset_flow_reconciles',sum(oldpatterns.values())==3857 and oldpatterns['complete']==result['n'] and 3857-result['n']==audit['old28_cases_removed_by_routine_n'])
    check('each_exclusion_reason_reconciles_overlapping_patterns',all(sum(v for pattern,v in patterns.items() if name in pattern.split('|'))==count for name,count in audit['exclusion_reason_counts'].items()))
    check('old28_patterns_only_routine_reasons',all(set(p.split('|'))<={'complete','BMI','HbA1c','eGFR','UACR'} for p in oldpatterns))
    check('routine_complete_then_protein_flow',audit['routine_complete_n']-sum(v for p,v in patterns.items() if p=='protein')==result['n'])
    check('source_flag_totals_common',all(sum(c.values())==result['n'] for c in audit['source_flag_counts_in_common'].values()))
    check('source_flag_totals_eligible',all(sum(audit[k].values())==3899 for k in ('glyco_method_counts','bmi_weight_flag_counts','bmi_height_flag_counts')))
    check('suspect_urine_and_fill_flags_reconcile',audit['urine_creatinine_suspect_replacement_n']==audit['urine_creatinine_below10_n'] and audit['common_urine_albumin_point4_n']<=audit['urine_albumin_replacement_point4_n'])
    flagcounts=audit['source_flag_counts_in_common']
    strict_exclusions=[result['n']-flagcounts['BMPWTFLG'].get('0.0',0),result['n']-flagcounts['BMPHTFLG'].get('0.0',0),flagcounts['GHPMETH'].get('4.0',0),audit['common_urine_albumin_point4_n']]
    strict_removed=result['n']-audit['strict_metadata_domain_n']
    check('strict_exclusion_overlap_union_bounds',max(strict_exclusions)<=strict_removed<=sum(strict_exclusions),{'marginal_nonzero_weightflag_heightflag_method4_UBPfill':strict_exclusions,'union_removed':strict_removed})
    check('calibration_intercept_slope_defined_only_counts_explicit',all(result['primary_intervals'][name+'|'+key]['valid']==498 and result['primary_intervals'][name+'|'+key]['undefined']==2 and result['primary_intervals'][name+'|'+key]['defined_draws_only'] is True for name in PRIMARY for key in ('calibration_intercept','calibration_slope')))
    check('all_eleven_routine_features_ordered',result['routine_features']==['age_linear','age_quadratic','male','current','former','sbp','bmi','hba1c','egfr_low','egfr_high','log_uacr'])
    check('no_duplicate_creatinine_or_diabetes_predictor',all('creatinine' not in name and 'diabetes' not in name for name in result['routine_features']))
    check('all_six_prespecified_arms',set(result['point'])==ARMS)
    for name,m in result['point'].items():metric_checks('main_'+name,m,result['n'],result['events8'])
    check('main_observed_rates_identical',len({m['observed'] for m in result['point'].values()})==1)
    check('main_weight_sums_identical',len({m['weight_sum'] for m in result['point'].values()})==1)
    for key in ('auc','brier'):
        check('primary_delta_'+key+'_arithmetic',same(result['primary']['delta_'+key],result['point']['routine_all3'][key]-result['point']['routine'][key],1e-15))
    check('primary_contrast_and_horizon_locked',result['primary_contrast']=='routine_all3_minus_routine' and result['horizon_years']==8 and result['weight']=='WTCYPEX6' and result['plan_commit']=='7800458')

    domains=Counter(f['analysis_domain'] for f in fits)
    validpart=[r for r in result['partition_sensitivity'] if r['status']=='computed']
    strict=result['strict_metadata_sensitivity']
    expected={'main':30,**{'partition_'+str(r['partition']):10 for r in validpart}}
    if strict['status']=='computed':expected['strict_metadata']=10
    check('exported_fold_parameter_domains_complete',dict(domains)==expected,dict(domains))
    check('unique_domain_model_outer_fold',len({(f['analysis_domain'],f['model'],f['outer_fold']) for f in fits})==len(fits))
    check('whole_PSU_outer_overlap_zero',all(f['PSU_overlap_n']==0 for f in fits))
    check('positive_training_event_nonevent_support',all(0<f['training_positive_weight_events']<f['training_positive_weight_n'] for f in fits))
    check('normalized_training_weights_equal_raw_n',all(same(f['normalized_weight_sum'],f['training_n'],1e-8) for f in fits))
    check('convergence_status_gradient_consistent',all((f['convergence_status']=='gradient_1e-6' and f['gradient_max']<=1e-6) or (f['convergence_status']=='tiny_step_gradient_1e-4' and f['gradient_max']<=1e-4) for f in fits))
    check('all_penalties_on_locked_grid',all(f['lambda']==f['selected_penalty'] and f['lambda'] in (1.,10.,100.) and set(map(float,f['inner_brier_by_penalty']))=={1.,10.,100.} for f in fits))
    check('training_inner_Brier_selection_and_tie_rule',all(f['lambda']==max(float(p) for p,s in f['inner_brier_by_penalty'].items() if s<=min(f['inner_brier_by_penalty'].values())+1e-12) for f in fits))
    dims={'simple':6,'simple_all3':9,'routine':11,'routine_all3':14,'rank_all3':14,'clip_all3':14}
    check('feature_dimensions_and_scales',all(len(f['training_mean'])==len(f['training_scale'])==len(f['feature_names'])==dims[f['model']] and f['core_feature_count']==(6 if f['model'].startswith('simple') else 11) and f['core_feature_indices']==list(range(f['core_feature_count'])) and len(f['core_coefficients_standardized'])==f['core_feature_count'] and len(f['protein_coefficients_standardized'])==dims[f['model']]-f['core_feature_count'] and all(math.isfinite(v) and v>0 for v in f['training_scale']) for f in fits))
    check('constant_features_scale_one',all(all(not c or s==1 for c,s in zip(f['constant_features'],f['training_scale'])) for f in fits))
    check('no_empirical_rank_support_values_exported',all(all('values' not in t and 'cdf' not in t for t in f['protein_transform_parameters']) for f in fits))
    for domain in expected:
        fs=[f for f in fits if f['analysis_domain']==domain]
        names=ARMS if domain=='main' else PRIMARY
        n=strict['n'] if domain=='strict_metadata' else result['n']
        events=strict['events8'] if domain=='strict_metadata' else result['events8']
        check(domain+'_complete_paired_arm_folds',{(f['model'],f['outer_fold']) for f in fs}=={(m,k) for m in names for k in range(5)})
        check(domain+'_train_scoring_population_totals',all(f['training_n']+f['scoring_n']==n and f['training_events']+f['scoring_events']==events for f in fs))
        baseline={f['outer_fold']:tuple(f[k] for k in ('training_n','scoring_n','training_events','scoring_events','training_psu_n','scoring_psu_n')) for f in fs if f['model']==('simple' if domain=='main' else 'routine')}
        check(domain+'_identical_partitions_all_arms',all(tuple(f[k] for k in ('training_n','scoring_n','training_events','scoring_events','training_psu_n','scoring_psu_n'))==baseline[f['outer_fold']] for f in fs))
        check(domain+'_scoring_folds_cover_population',all(sum(f['scoring_n'] for f in fs if f['model']==name)==n and sum(f['scoring_events'] for f in fs if f['model']==name)==events for name in names))
    for name in ARMS:
        modelbins=[b for b in bins if b['model']==name and b['analysis_domain']=='main'];total=sum(b['weight_sum'] for b in modelbins)
        check(name+'_calibration_bins_population',sum(b['n'] for b in modelbins)==result['n'] and sum(b['events'] for b in modelbins)==result['events8'])
        check(name+'_calibration_bin_weighted_rates',same(sum(b['weight_sum']*b['observed'] for b in modelbins)/total,result['point'][name]['observed']) and same(sum(b['weight_sum']*b['predicted'] for b in modelbins)/total,result['point'][name]['predicted']) and same(total,result['point'][name]['weight_sum'],1e-6))
    for name in PRIMARY:check(name+'_bootstrap_penalty_counts',sum(result['bootstrap_penalty_selection_counts'][name].values())==5*(500-failures))
    partitions=result['partition_sensitivity']
    check('all_five_partitions_retained',len(partitions)==5 and sorted(r['partition'] for r in partitions)==list(range(5)))
    check('partition_failures_explicit',result['partition_failures']==5-len(validpart))
    check('partition_sensitivity_separate_from_CI',result['partition_sensitivity_pooled_into_CI'] is False)
    for key in ('auc','brier'):
        contrast='delta_'+key
        check('partition_'+contrast+'_arithmetic',all(same(r[contrast],r['routine_all3|'+key]-r['routine|'+key],1e-15) for r in validpart))
        check('partition_'+contrast+'_range',same(result['partition_sensitivity_range'][contrast]['minimum'],min(r[contrast] for r in validpart)) and same(result['partition_sensitivity_range'][contrast]['maximum'],max(r[contrast] for r in validpart)))
    check('strict_domain_nested_and_conditional',strict['n']==audit['strict_metadata_domain_n'] and 0<strict['n']<=result['n'] and strict['interpretation']=='different_conditional_domain_within_domain_pair_only')
    if strict['status']=='computed':
        check('strict_primary_pair_same_domain',set(strict['point'])==PRIMARY)
        for name,m in strict['point'].items():metric_checks('strict_'+name,m,strict['n'],strict['events8'])
        for key in ('auc','brier'):check('strict_delta_'+key+'_arithmetic',same(strict['primary']['delta_'+key],strict['point']['routine_all3'][key]-strict['point']['routine'][key],1e-15))
        check('strict_no_bootstrap_CI_claimed','intervals' not in strict and 'primary_intervals' not in strict)
    check('research_limitations_explicit',all(result[k] is False for k in ('assay_and_creatinine_calibration_uncertainty_included','reported_protein_QC_LOQ_validated','independent_external_validation','clinical_use_ready','participant_records_exported')) and result['same_common_domain_all_primary_and_simple_comparators'] is True and result['no_single_protein_winner_search'] is True)
    limits=[
        'No participant records, raw source archives or individual predictions were accessed, and no model was rerun. Point Brier/AUC cannot be independently re-estimated from these aggregate exports.',
        'All 500 paired contrast rows were parsed independently and percentile endpoints recomputed by sorting and linear interpolation with Python standard library only.',
        'Calibration intercept/slope intervals for both primary arms have only 498 defined draws and 2 undefined; their endpoints were matched to stdout but cannot be independently recomputed from the primary-only CSV. They are conditional defined-draw intervals, not nominal full-500 coverage.',
        'Source fingerprints, SAS hashes and reader aggregate counts were matched to the successful job manifest; source data bytes were not downloaded locally. Additional direct official lab.sas metadata GET was denied by proxy (Tunnel connection failed: 403 Forbidden), so no independent local source-byte verification is claimed.',
        'Exported fit parameters cover main, partition and strict point refits. Per-bootstrap fit parameters are not exported; model convergence guards were reviewed and exported failure/penalty counts checked.',
        'Primary intervals condition on the fixed original PSU partition, complete profiles, reported assay values and specified model forms. Five additional partitions are exploratory sensitivities and are not pooled into those intervals.',
        'NHANES III is reused; this is internal research, with remaining historical assay/QC/LOQ and creatinine calibration uncertainty. Neither independent external validation nor clinical readiness is established.'
    ]
    report={'audit':'independent_routine29_aggregates_provenance_and_refit_percentiles','run_id':RUN,'job_id':JOB,'checks':checks,'details':details,'passed':sum(checks.values()),'total':len(checks),'failed_checks':errors,'scientific_python_files_checked':len(scientific),'result_sha256':digest(ROOT/'RESULTS_SUMMARY.json'),'primary_bootstrap_csv_sha256':digest(ROOT/'primary_bootstrap.csv'),'independent_primary_limits':independent_limits,'common_n':result['n'],'events8':result['events8'],'limits':limits}
    (OUT/'verification.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    lines=['# Independent routine29 aggregate and percentile audit','',f"{report['passed']}/{report['total']} checks passed; {len(scientific)} scientific Python files checked.",'',f"Common n={result['n']}, events8={result['events8']}; primary independent percentiles: {json.dumps(independent_limits)}.",'',*['- '+s for s in limits]]
    if errors:lines+=['','Failed checks: '+', '.join(errors)]
    (OUT/'verification.md').write_text('\n'.join(lines)+'\n')
    print(json.dumps({k:report[k] for k in ('passed','total','failed_checks','scientific_python_files_checked','independent_primary_limits','result_sha256')},indent=2))
    return report

if __name__=='__main__':run()
