"""Independent stdlib audit of freeze30 aggregates and synthetic scoring only.

No participant data, raw downloads, training or model refits. The reference
feature equations, standardization and expit below use Python's math module.
A subprocess calls the reviewed adapter only on invented numeric fixtures.
"""
from pathlib import Path
import argparse, hashlib, importlib.util, json, math, re, subprocess, sys

REPO=Path(__file__).resolve().parents[1];PARENT=REPO/'research/mortality';ROOT=PARENT/'freeze30'
EVIDENCE=Path('/tmp');OUT=Path('/tmp');LOG=EVIDENCE/'freeze30_independent_job.log';RUN=40903512074;JOB=123551907162
SCIENCE_COMMIT='16e2eb8fe0abd51b351a6c3775a7180ebdd25354'
ARMS=('routine','routine_all3');checks={};details={};errors=[]

def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def check(name,value,detail=None):
    checks[name]=bool(value)
    if detail is not None:details[name]=detail
    if not value:errors.append(name)

def same(a,b,tol=1e-12):
    return a is not None and b is not None and math.isfinite(float(a)) and math.isfinite(float(b)) and abs(float(a)-float(b))<=tol

def blocks(text,label):
    begin='FREEZE30_'+label+'_BEGIN';end='FREEZE30_'+label+'_END';result=[];active=None
    for line in text.splitlines():
        if line==begin:
            if active is not None:raise ValueError('Nested marker')
            active=[]
        elif line==end:
            if active is None:raise ValueError('Unopened marker')
            result.append('\n'.join(active)+'\n');active=None
        elif active is not None:active.append(line)
    if active is not None:raise ValueError('Unterminated marker')
    return result

def converged(row):
    return ((row['convergence_status']=='gradient_1e-6' and 0<=row['gradient_max']<=1e-6)
        or (row['convergence_status']=='tiny_step_gradient_1e-4' and 0<=row['gradient_max']<=1e-4))

def allkeys(value):
    if isinstance(value,dict):
        for key,item in value.items():yield key;yield from allkeys(item)
    elif isinstance(value,list):
        for item in value:yield from allkeys(item)

def finite(value):return type(value) in (int,float) and math.isfinite(value)

def fixture(**changes):
    row=dict(age_years=50,sex='male',smoking='current',sbp_mmhg=150.,sbp_measurement_count=3,
        bmi_kg_m2=32.,bmi_weight_flag=0,bmi_height_flag=0,hba1c_percent=6.5,hba1c_method=1,
        creatinine_historical_mg_dl=(1.+.184)/.960,urine_creatinine_mg_dl=100.,urine_albumin_ug_ml=10.,
        cystatin_c_reported_mg_l=.5,beta2_microglobulin_reported_mg_l=2.,beta_trace_protein_reported_mg_l=4.)
    row.update(changes);return row

def independent_features(row):
    female=row['sex']=='female';age=(row['age_years']-60.)/10.
    serum=.960*row['creatinine_historical_mg_dl']-.184
    ratio=serum/(.7 if female else .9)
    egfr=142.*min(ratio,1.)**(-.241 if female else -.302)*max(ratio,1.)**(-1.2)*.9938**row['age_years']*(1.012 if female else 1.)
    uacr=100.*row['urine_albumin_ug_ml']/row['urine_creatinine_mg_dl']
    core=[age,age*age,float(not female),float(row['smoking']=='current'),float(row['smoking']=='former'),
        (row['sbp_mmhg']-130.)/20.,(row['bmi_kg_m2']-27.)/5.,row['hba1c_percent'],
        min(egfr,60.)/15.,max(egfr-60.,0.)/30.,math.log2(uacr)]
    protein=[math.log2(row[name]) for name in ('cystatin_c_reported_mg_l','beta2_microglobulin_reported_mg_l','beta_trace_protein_reported_mg_l')]
    return core,protein,egfr,uacr

def expit(value):
    if value>=0:return 1./(1.+math.exp(-value))
    e=math.exp(value);return e/(1.+e)

def run():
    raw=LOG.read_text();plain='\n'.join(re.sub(r'^\S+ ','',line) for line in raw.splitlines())
    payload={}
    for label,name in [('ARTIFACT','FROZEN_MODELS.json'),('MANIFEST','FREEZE_MANIFEST.json')]:
        found=blocks(plain,label)
        check(label+'_exactly_two_successful_run_blocks',len(found)==2)
        if len(found)!=2:raise ValueError('Expected two successful run outputs')
        check(label+'_two_run_bytes_identical',found[0]==found[1])
        check(label+'_local_bytes_equal_both_stdout_runs',(ROOT/name).read_bytes()==found[0].encode())
        payload[label]=json.loads(found[0])
    artifact=payload['ARTIFACT'];manifest=payload['MANIFEST']
    check('artifact_hash_matches_manifest',digest(ROOT/'FROZEN_MODELS.json')==manifest['artifact_sha256'])
    check('artifact_expected_scientific_run_hash',digest(ROOT/'FROZEN_MODELS.json')=='70d29bf6e60d0967a11534634ebd1769a8b513b97ecb1fe8c66c6a5392534d43')
    check('manifest_expected_successful_run_hash',digest(ROOT/'FREEZE_MANIFEST.json')=='f33df537280ead416406139412c4155f05ba8215e37408873f8442df7c4ab488')
    check('artifact_and_manifest_schema',artifact['schema']=='freeze30-model/1' and manifest['schema']=='freeze30-manifest/1' and manifest['artifact_file']=='FROZEN_MODELS.json')
    lock=json.loads((ROOT/'LOCK.json').read_text());contract=json.loads((ROOT/'CONTRACT.json').read_text())
    for field,name in [('contract_sha256','CONTRACT.json'),('plan_sha256','FREEZE_PLAN.md'),('lock_sha256','LOCK.json')]:
        check(field+'_artifact_manifest_current_file',artifact[field]==manifest[field]==digest(ROOT/name))
    check('all_locked_metadata_hashes_match',all(digest(PARENT/name if '/' in name else ROOT/name)==expected for name,expected in lock['metadata_sha256'].items()))
    check('plan_ancestor_of_scientific_code',subprocess.run(['git','merge-base','--is-ancestor','5cbf90e',SCIENCE_COMMIT],cwd=REPO,capture_output=True).returncode==0)
    check('plan_commit_fixed_in_artifact',artifact['plan_commit']=='5cbf90e')
    scientific=artifact['scientific_python_sha256'];old=json.loads((PARENT/'routine29/RESULTS_SUMMARY.json').read_text())
    check('artifact_manifest_scientific_hash_registry_equal',scientific==manifest['scientific_python_sha256'])
    check('all_84_scientific_hashes_match_current_files',len(scientific)==84 and all(digest(PARENT/name)==h for name,h in scientific.items()))
    check('every_current_scientific_py_is_recorded',set(scientific)=={str(p.relative_to(PARENT)) for p in PARENT.rglob('*.py')})
    check('all_80_historical29_hashes_preserved',len(old['scientific_python_sha256'])==80 and all(scientific[name]==h for name,h in old['scientific_python_sha256'].items()))
    check('only_four_new_scientific_files',set(scientific)-set(old['scientific_python_sha256'])=={'freeze30/__init__.py','freeze30/adapter.py','freeze30/freeze.py','freeze30/qualification.py'})
    changed=subprocess.check_output(['git','diff','--name-status',lock['base_commit'],'--','research/mortality','tests'],cwd=REPO,text=True).splitlines()
    check('old_mortality_and_test_files_unchanged_from_base',all(line.startswith('A\t') and line.split('\t')[1].startswith(('research/mortality/freeze30/','tests/test_freeze30')) for line in changed))
    historical=json.loads((PARENT/'temporal26/LOCK.json').read_text())['models']
    check('all_historical_frozen_model_hashes_preserved',all(digest(PARENT/name)==h for name,h in historical.items()))

    jobs=json.loads((EVIDENCE/'freeze30_jobs.json').read_text())['jobs']
    check('three_CI_jobs_completed_success',len(jobs)==3 and all(j['status']=='completed' and j['conclusion']=='success' for j in jobs))
    for version,file in [('3.11',EVIDENCE/'freeze30_synthetic311.log'),('3.13',EVIDENCE/'freeze30_synthetic313.log')]:
        text=Path(file).read_text();check('CI_'+version+'_60_tests_pass',re.search(r'Ran 60 tests in',text) is not None and re.search(r'(?m)^\S+ OK$',text) is not None)
    check('successful_job_contains_two_byte_comparison_commands','cmp /tmp/freeze30_first/FROZEN_MODELS.json /tmp/freeze30_replay/FROZEN_MODELS.json' in plain and 'cmp /tmp/freeze30_first/FREEZE_MANIFEST.json /tmp/freeze30_replay/FREEZE_MANIFEST.json' in plain)

    for field in ('source_manifest_sha256','source_content_sha256','source_fingerprints','source_policy','training'):
        check(field+'_artifact_manifest_equal',artifact[field]==manifest[field])
    fingerprints=artifact['source_fingerprints'];prior=old['source_audit']['routine_source_fingerprints'];policy=artifact['source_policy'];audit=manifest['source_audit']
    check('14_source_file_names_unchanged',set(fingerprints)==set(prior) and len(fingerprints)==14)
    check('all_13_nonindex_source_hashes_unchanged',all(fingerprints[name]==h for name,h in prior.items() if name!='index.htm'))
    check('only_index_source_hash_changed',[name for name,h in fingerprints.items() if h!=prior[name]]==['index.htm'])
    check('index_change_documented_exactly',policy['pinned_non_index_files']==13 and policy['index_is_dynamic_catalog_not_scientific_table'] is True and policy['prior_index_sha256']==prior['index.htm'] and policy['actual_index_sha256']==fingerprints['index.htm'] and policy['index_changed'] is True and policy['every_non_index_source_sha_unchanged'] is True and policy['every_resolved_release_url_unchanged'] is True)
    check('source_audit_fingerprint_registry_matches',audit['routine_source_fingerprints']==fingerprints)
    stable=hashlib.sha256(json.dumps(fingerprints,sort_keys=True,separators=(',',':')).encode()).hexdigest()
    check('stable_fingerprint_content_hash_recomputed',artifact['source_content_sha256']==stable==audit['routine_source_content_sha256'])
    check('retrieval_manifest_hash_consistent_with_source_audit',artifact['source_manifest_sha256']==audit['routine_source_manifest_sha256'])
    for kind,n,fields,cells,end,width in [('lab',29314,6,175884,1970,1977),('exam',31311,4,125244,1533,6233)]:
        reader=audit['routine_readers'][kind];previous=old['source_audit']['routine_readers'][kind]
        check(kind+'_reader_aggregates_exactly_preserved',reader==previous)
        check(kind+'_independent_reader_cells',reader['exact_match'] is True and reader['rows']==n and reader['fields']==fields and reader['cell_comparisons']==cells)
        check(kind+'_source_and_SAS_hashes_agree',reader['source_sha256']==fingerprints[kind+'.dat'] and reader['sas_sha256']==fingerprints[kind+'.sas'])
        check(kind+'_width_and_scaling_contract',reader['SAS_MISSOVER'] is True and reader['selected_fields_end']==end and reader['nominal_content_width']==width and reader['FORMAT_scaling_applied'] is False and sum(reader['source_record_width_counts'].values())==n)
    check('common_domain_audit_preserved',all(audit[k]==old['source_audit'][k] for k in ('eligible_protein_weight_n','old28_common_n','routine_common_n','old28_cases_removed_by_routine_n','exclusion_reason_patterns','source_flag_counts_in_common','missing_routine')))
    training=artifact['training']
    check('training_domain_n_events_fixed',training['n']==lock['expected_training_n']==old['n']==3730 and training['events8']==lock['expected_events8']==old['events8']==563)
    check('training_weights_design_preserved',training['positive_weight_n']==3730 and training['weight']=='WTCYPEX6' and training['strata_n']==49 and training['psu_n']==98 and same(training['effective_n'],old['effective_n']) and same(training['weight_sum'],old['point']['routine']['weight_sum'],1e-6))
    check('final_inner_algorithm_locked',training['inner_folds']==lock['final_inner_folds']==3 and training['inner_salt']==lock['inner_salt']=='freeze30:final-inner:v1' and training['weight_normalization']=='sum_to_raw_training_n')
    check('target_matches_contract',artifact['target']==contract['target'] and artifact['target']['horizon_months']==96 and artifact['target']['time_origin']=='first_exam' and artifact['target']['endpoint']=='all_cause_death')
    check('historical_only_mapping_and_no_external_approval',artifact['mapping_id']==contract['mapping_id']=='NHANESIII_FIRST_EXAM_REPORTED_V1' and contract['external_mappings_approved']==[])
    check('model_family_and_training_scope_explicit',artifact['model_family']=='eight_year_survey_weighted_ridge_logistic' and artifact['fit_scope']=='full_training_final_fit_not_new_validation' and artifact['selection_scope']=='training_only_not_validation')
    check('all_literal_research_flags_false',all(artifact[k] is False for k in ('clinical_use_ready','independent_external_validation','participant_records_exported')) and all(manifest[k] is False for k in ('clinical_use_ready','training_performance_exported','participant_records_exported')))
    forbidden={'records','SEQNs','SEQN','participant_id','participants','training_predictions','predictions','outcomes','OOF','oof_predictions','training_auc','training_brier'}
    check('no_record_ID_outcome_or_prediction_keys_exported',not (set(allkeys(artifact))|set(allkeys(manifest)))&forbidden)
    check('exactly_two_prespecified_models',set(artifact['models'])==set(ARMS))
    names=['age_centered','age_centered_squared','male','smoking_current','smoking_former','systolic_bp','bmi','hba1c','egfr_low','egfr_high','log2_uacr']
    for arm in ARMS:
        row=artifact['models'][arm];prep=row['preprocessing'];fit=row['fit'];dim=11 if arm=='routine' else 14
        check(arm+'_feature_names_order',row['feature_names']==names+([] if arm=='routine' else ['log2_cystatin_c','log2_beta2_microglobulin','log2_beta_trace_protein']))
        check(arm+'_feature_dimension_and_finiteness',all(isinstance(x,list) and len(x)==dim and all(finite(v) for v in x) for x in (prep['mean'],prep['scale'],fit['coefficients'])) and all(v>0 for v in prep['scale']) and finite(fit['intercept']))
        check(arm+'_constant_feature_flags',len(prep['constant'])==dim and all(type(v) is bool for v in prep['constant']) and all(not c or s==1 for c,s in zip(prep['constant'],prep['scale'])))
        check(arm+'_exact_frozen_transform',prep['model_name']==arm and prep['core_indices']==list(range(11)) and prep['protein_transforms']==([] if arm=='routine' else [{'column':k,'kind':'identity_log2'} for k in range(3)]))
        scores=row['inner_brier_by_penalty'];check(arm+'_fixed_penalty_grid',set(map(float,scores))=={1.,10.,100.} and all(finite(v) and 0<=v<=1 for v in scores.values()))
        minimum=min(scores.values());chosen=max(float(p) for p,score in scores.items() if score<=minimum+1e-12)
        check(arm+'_lambda_independent_selection_arithmetic',chosen==row['selected_penalty']==fit['penalty'],{'selected':chosen,'scores_are_training_diagnostics_only':True})
        check(arm+'_final_gradient_convergence',converged(fit) and finite(fit['objective']) and fit['objective']>0 and type(fit['iterations']) is int and fit['iterations']>0)
        check(arm+'_final_normalized_weight_sum',same(fit['normalized_weight_sum'],3730))
        fs=row['inner_folds_aggregate'];ds=row['inner_fit_diagnostics']
        check(arm+'_three_unique_inner_folds',len(fs)==3 and sorted(f['fold'] for f in fs)==[0,1,2])
        check(arm+'_inner_counts_cover_full_population',sum(f['scoring_n'] for f in fs)==3730 and sum(f['scoring_events'] for f in fs)==563 and sum(f['scoring_psu_n'] for f in fs)==98)
        check(arm+'_inner_training_scoring_partition',all(f['training_n']+f['scoring_n']==3730 and f['training_events']+f['scoring_events']==563 and f['training_psu_n']+f['scoring_psu_n']==98 and f['PSU_overlap_n']==0 for f in fs))
        check(arm+'_inner_positive_training_event_nonevent_support',all(0<f['training_positive_weight_events']<f['training_positive_weight_n']==f['training_n'] and f['scoring_positive_weight_n']==f['scoring_n'] for f in fs))
        check(arm+'_nine_unique_inner_optimizations',len(ds)==9 and {(d['fold'],d['penalty']) for d in ds}=={(f,p) for f in range(3) for p in (1.,10.,100.)})
        check(arm+'_inner_gradient_convergence',all(converged(d) and type(d['iterations']) is int and d['iterations']>0 for d in ds))
        check(arm+'_inner_weight_normalization',all(same(d['normalized_weight_sum'],next(f['training_n'] for f in fs if f['fold']==d['fold'])) for d in ds))
        replay=manifest['numerical_replay'][arm];check(arm+'_runner_training_adapter_replay_aggregate',replay['compared_n']==3730 and replay['maximum_absolute_probability_difference']==0.)
    check('paired_models_identical_inner_participants_folds',artifact['models']['routine']['inner_folds_aggregate']==artifact['models']['routine_all3']['inner_folds_aggregate'])
    check('paired_final_routine_scaling_identical',all(artifact['models']['routine']['preprocessing'][k]==artifact['models']['routine_all3']['preprocessing'][k][:11] for k in ('mean','scale','constant')))

    # Invented fixtures only; do not import scientific feature engineering for
    # the reference calculations or run any fit or optimizer.
    rows=[fixture(),fixture(age_years=40,smoking='never',urine_creatinine_mg_dl=10.),
        fixture(age_years=75,sex='female',smoking='former',creatinine_historical_mg_dl=2.,urine_albumin_ug_ml=.4),
        fixture(age_years=60,sex='female',smoking='never',sbp_mmhg=130.,bmi_kg_m2=27.,hba1c_percent=5.4,hba1c_method=4,bmi_weight_flag=3,bmi_height_flag=2),
        fixture(age_years=70,sex='male',creatinine_historical_mg_dl=5.,urine_albumin_ug_ml=500.,cystatin_c_reported_mg_l=2.,beta2_microglobulin_reported_mg_l=4.,beta_trace_protein_reported_mg_l=.25),
        fixture(age_years=50,sex='female',creatinine_historical_mg_dl=.4,urine_albumin_ug_ml=.4,urine_creatinine_mg_dl=100.)]
    document={'schema':contract['input_schema'],'contract_sha256':digest(ROOT/'CONTRACT.json'),'mapping_id':contract['mapping_id'],'units':contract['units'],'records':rows}
    adapter_code="""import json,sys
from research.mortality.freeze30 import adapter
doc=json.load(sys.stdin)
artifact=adapter.load_artifact('research/mortality/freeze30/FROZEN_MODELS.json','research/mortality/freeze30/FREEZE_MANIFEST.json')
out=adapter.score(doc,artifact)
print(json.dumps({key:value.tolist() for key,value in out.items()},allow_nan=False))
"""
    executed=subprocess.run([sys.executable,'-c',adapter_code],cwd=REPO,input=json.dumps(document),text=True,capture_output=True,check=True)
    actual=json.loads(executed.stdout);synthetic_difference={}
    check('independent_known_eGFR_equation_fixture',same(independent_features(rows[0])[2],91.69147860983348,1e-12))
    check('independent_urine_albumin_fill_fixture',same(independent_features(rows[-1])[3],.4,1e-15))
    for arm in ARMS:
        row=artifact['models'][arm];prep=row['preprocessing'];fit=row['fit'];expected=[];raw_expected=[]
        raw_beta=[b/s for b,s in zip(fit['coefficients'],prep['scale'])]
        raw_intercept=fit['intercept']-sum(m*b for m,b in zip(prep['mean'],raw_beta))
        for invented in rows:
            core,protein,_,_=independent_features(invented);values=core if arm=='routine' else core+protein
            eta=fit['intercept']+sum((x-m)/s*b for x,m,s,b in zip(values,prep['mean'],prep['scale'],fit['coefficients']))
            expected.append(expit(eta));raw_expected.append(expit(raw_intercept+sum(x*b for x,b in zip(values,raw_beta))))
        diff=max(abs(a-b) for a,b in zip(actual[arm],expected));synthetic_difference[arm]=diff
        check(arm+'_six_synthetic_scores_match_independent_stdlib_math',len(actual[arm])==6 and diff<=1e-12,{'synthetic_n':6,'maximum_absolute_difference':diff})
        check(arm+'_raw_coefficient_reconstruction_equivalence',all(same(a,b,1e-12) for a,b in zip(raw_expected,expected)))

    # Pure metadata replay, separate from and incapable of participant scoring.
    spec=importlib.util.spec_from_file_location('independent_review_qualification30',ROOT/'qualification.py');q=importlib.util.module_from_spec(spec);spec.loader.exec_module(q)
    candidates=json.loads((ROOT/'CANDIDATE_MANIFEST.json').read_text());qualification=json.loads((ROOT/'QUALIFICATION_RESULTS.json').read_text())
    replay=[dict(candidate_id=c['candidate_id'],**q.audit(c['qualification_metadata'])) for c in candidates['candidates']]
    check('qualification_12_result_metadata_replay_exact',replay==qualification['candidates'] and len(replay)==qualification['candidate_count']==12)
    check('qualification_every_candidate_blocked',all(r['status']=='blocked_before_scoring' and r['frozen_scoring_allowed'] is False for r in replay))
    check('qualification_no_external_candidates_comparisons_or_files',qualification['eligible_external_candidates']==0 and qualification['new_external_outcome_comparisons']==0 and qualification['new_external_participant_files_retrieved']==0 and qualification['external_mappings_approved']==[] and qualification['clinical_use_ready'] is False)
    check('external_mapping_registry_empty',q.QUALIFIED_EXTERNAL_MAPPINGS=={})
    check('qualification_does_not_authenticate_access_independence_or_chronology',all(all(r[k] is False for k in ('legal_authorization_verified','source_authority_verified','independence_verified','chronology_verified','person_records_processed','network_or_file_access_performed')) for r in replay))
    ids={r['candidate_id']:r for r in replay}
    check('HRS_AI_policy_blocker_preserved','HRS_person_records_prohibited_in_AI_or_unattended_cloud' in ids['HRS2016VBS']['access_blockers'])
    check('NHANESIII_first_and_repeat_reuse_blocked',all('NHANESIII_first_second_or_renamed_release_already_used' in ids[name]['overlap_blockers'] for name in ('NHANESIII_FIRST_EXAM','NHANESIII_SECOND_EXAM')))
    check('SSCARD_old_calendar_cohort_reuse_blocked','NHANES1999_2018_programme_cycle_already_used' in ids['NHANES1999_2004_SSCARD']['overlap_blockers'])
    sscard=json.loads((ROOT/'SSCARD_SOURCE_METADATA.json').read_text())
    check('SSCARD_history_and_missing_BTP_explicit',sscard['previously_used_participant_cohort'] is True and sscard['prior_beta2_loaded'] is True and sscard['independent_new_cohort'] is False and sscard['full_three_protein_source_gate_pass'] is False and sscard['beta_trace_protein_present'] is False)
    limits=[
        'No participant records, raw archive requests or model refits were performed. Whole-PSU membership and training fit equations were reviewed in code; aggregate counts alone cannot independently reconstruct participant-level splits, gradients or training Brier.',
        'The two artifacts and manifests were compared byte-for-byte against both outputs of an independently fetched successful CI job. This is same-environment deterministic replay, not a cross-platform guarantee.',
        'The reference scoring equations and raw-coefficient reconstruction use stdlib math on six explicitly invented fixtures. The adapter was called only on those synthetic fixtures; no training participants were scored by this reviewer.',
        'All13 scientific/documentation source hashes match the preserved29 hashes, and the recorded current source digest was independently recomputed. Source bytes were not downloaded locally. The new full retrieval-manifest bytes are not exported, so its SHA was cross-checked between artifact/source audit/manifest rather than independently rehashed from that unexported file.',
        'Zero maximum training replay differences are successful-run aggregate assertions, not local participant-level replay.',
        'Qualification results replay metadata assertions only. Empty external mapping registry, zero qualified candidates, historical assay limits and lawful-access requirements remain; these artifacts are not independent external or clinical validation.'
    ]
    report={'audit':'independent_freeze30_aggregates_hashes_and_stdlib_synthetic_scoring','scientific_commit':SCIENCE_COMMIT,'run_id':RUN,'job_id':JOB,'checks':checks,'details':details,'passed':sum(checks.values()),'total':len(checks),'failed_checks':errors,'scientific_python_files_checked':len(scientific),'historical_scientific_hashes_preserved':80,'artifact_sha256':digest(ROOT/'FROZEN_MODELS.json'),'manifest_sha256':digest(ROOT/'FREEZE_MANIFEST.json'),'synthetic_fixture_n':6,'synthetic_max_abs_difference':synthetic_difference,'limits':limits}
    (OUT/'INDEPENDENT_VERIFICATION.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    lines=['# Independent freeze30 aggregate and synthetic-scoring audit','',f"{report['passed']}/{report['total']} checks passed;84 scientific files checked and80 historical hashes preserved.",'',f"Artifact SHA256 `{report['artifact_sha256']}`.",f"Manifest SHA256 `{report['manifest_sha256']}`.",'',f"Synthetic reference maximum absolute differences: {json.dumps(synthetic_difference)}.",'',*['- '+text for text in limits]]
    if errors:lines+=['','Failed checks: '+', '.join(errors)]
    (OUT/'INDEPENDENT_VERIFICATION.md').write_text('\n'.join(lines)+'\n')
    print(json.dumps({k:report[k] for k in ('passed','total','failed_checks','artifact_sha256','manifest_sha256','synthetic_max_abs_difference')},indent=2))
    return report

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--evidence-dir',type=Path,required=True,help='Directory containing four saved authorized job-log/metadata files')
    parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args();EVIDENCE=args.evidence_dir;OUT=args.out
    LOG=EVIDENCE/'freeze30_independent_job.log'
    OUT.mkdir(parents=True,exist_ok=False)
    result=run()
    raise SystemExit(1 if result['failed_checks'] else 0)
