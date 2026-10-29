"""Independent stdlib audit of B2M31 exported aggregates; never opens source records."""
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

CODE = 'f27c7b47157db7fbaf31dc998bbb5c55c1c4f7e8'
MERGE = 'a752d10a251f4bc0a381fe44d9d9c0a117d64f6e'


def digest(data):
    return hashlib.sha256(data).hexdigest()


def pairs(items):
    out = {}
    for key, value in items:
        if key in out:
            raise ValueError('Duplicate JSON key')
        out[key] = value
    return out


def strict(data):
    def reject(value):
        raise ValueError('Nonfinite JSON literal')
    out = json.loads(data, object_pairs_hook=pairs, parse_constant=reject)
    def finite(value):
        if isinstance(value,float) and not math.isfinite(value):
            raise ValueError('Nonfinite JSON numeric value')
        if isinstance(value,dict):
            for x in value.values():finite(x)
        elif isinstance(value,list):
            for x in value:finite(x)
    finite(out)
    return out


def safe_url(value):
    p = urllib.parse.urlsplit(value)
    return (p.scheme == 'https' and p.hostname in {'wwwn.cdc.gov', 'www.cdc.gov', 'ftp.cdc.gov'}
            and p.username is None and p.password is None and p.port in (None, 443)
            and bool(p.path) and not p.fragment)


def close(a, b):
    return math.isclose(a, b, rel_tol=1e-12, abs_tol=1e-9)


def audit(repo, log_path, ref_path):
    base = repo / 'research/mortality'
    folder = base / 'b2m31'
    result_bytes = (folder / 'RESULTS_SUMMARY.json').read_bytes()
    manifest_bytes = (folder / 'SOURCE_MANIFEST.json').read_bytes()
    results, manifest = strict(result_bytes), strict(manifest_bytes)
    lock_bytes = (folder / 'LOCK.json').read_bytes()
    lock = strict(lock_bytes)
    spec = strict((folder / 'SOURCE_SPEC.json').read_bytes())
    folds = strict((folder / 'FOLD_PARAMETERS.json').read_bytes())
    calibration = strict((folder / 'CALIBRATION_BINS.json').read_bytes())
    reference = strict(ref_path.read_bytes())
    raw_log = log_path.read_bytes()
    log = raw_log.decode('utf-8')
    checks = {}

    def require(name, condition):
        checks[name] = bool(condition)
        if not condition:
            raise AssertionError(name)

    def git_bytes(commit, path):
        return subprocess.check_output(['git', 'show', commit + ':' + path], cwd=repo)

    require('independent_pre_result_reference_hash', digest(ref_path.read_bytes())
            == 'd46459d8b5e43390a208153be53482e89b3c9be964a2efec9913c1a191aadc19')
    require('protocol_commit', manifest['plan_commit'] == results['plan_commit'] == '6774fbe')
    require('manifest_schema', manifest['schema'] == 'b2m31-retrieval/1')
    require('manifest_original_bytes_hash', digest(manifest_bytes) == results['source_manifest_export_sha256']
            == results['source_audit']['source_manifest_sha256'])
    require('lock_bound_to_result_and_manifest', digest(lock_bytes) == results['lock_sha256'] == manifest['lock_sha256'])
    require('plan_bound_to_result', digest((folder/'ANALYSIS_PLAN.md').read_bytes()) == results['plan_sha256'])
    require('twelve_locked_metadata_objects', len(lock['metadata_sha256']) == 12)
    metadata = {}
    for name, expected in lock['metadata_sha256'].items():
        p = base/name if '/' in name else folder/name
        rel = p.relative_to(repo).as_posix()
        observed = digest(p.read_bytes())
        require('locked_metadata:' + name, observed == expected == digest(git_bytes(CODE, rel)))
        metadata[name] = observed
    historical_code = {}
    for name, expected in lock['historical_scientific_python_sha256'].items():
        p = base/name
        require('historical_code:' + name, digest(p.read_bytes()) == expected)
        historical_code[name] = expected
    scientific = results['scientific_python_sha256']
    require('ninety_exported_scientific_files', len(scientific) == 90)
    for name, expected in scientific.items():
        p = base/name
        require('executed_science_hash:' + name, digest(p.read_bytes()) == expected
                == digest(git_bytes(CODE, p.relative_to(repo).as_posix())))
    require('historical_code_preserved_in_science_export',
            all(scientific.get(name) == expected for name, expected in historical_code.items()))

    files = manifest['files']
    rows = {r['file']: r for r in files}
    require('exact_fifty_unique_required_source_names', len(rows) == len(files) == 50
            and set(rows) == set(reference['required_files']))
    for name, row in rows.items():
        require('source_manifest_entry:' + name,
                Path(name).name == name and row['status'] == 'retrieved'
                and isinstance(row['bytes'], int) and not isinstance(row['bytes'], bool)
                and 0 < row['bytes'] <= 25_000_000
                and bool(re.fullmatch('[0-9a-f]{64}', row['sha256']))
                and safe_url(row['url']) and safe_url(row['resolved_url']))
    fingerprints = {name: row['sha256'] for name, row in rows.items()}
    content_hash = digest(json.dumps(fingerprints, sort_keys=True, separators=(',', ':')).encode())
    a = results['source_audit']
    require('fifty_exported_fingerprints_match_manifest', fingerprints == a['source_fingerprints']
            and a['files_verified'] == 50 and a['source_content_sha256'] == content_hash)
    data = {n: r for n, r in rows.items() if n.endswith(('.xpt', '.dat'))}
    require('exact_twentyfive_historical_data_hashes', len(data) == a['historical_data_hashes_verified'] == 25
            and {n:r['sha256'] for n,r in data.items()} == reference['historical_data_sha256'])
    sscard = rows['1999_SSCARD_A.xpt']
    expected_url = 'https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/1999/DataFiles/SSCARD_A.xpt'
    expected_urls = {}
    for cycle in spec['cycles']:
        for component, info in cycle['components'].items():
            expected_urls[str(cycle['year'])+'_'+component+'.xpt'] = info['historical_declared_data_url']
            expected_urls[str(cycle['year'])+'_'+component+'.htm'] = info['documentation_url']
        expected_urls[str(cycle['year'])+'_mortality.dat'] = cycle['mortality_url']
    expected_urls.update({'index.htm':spec['sscard_index_url'], '1999_SSCARD_A.htm':spec['sscard']['metadata_url'],
                          '1999_SSCARD_A.xpt':expected_url,'weighting.htm':spec['weighting_url'], 'variance.htm':spec['variance_url']})
    require('all_fifty_declared_urls_match_pinned_spec_or_runtime_sscard',
            {n:r['url'] for n,r in rows.items()} == expected_urls)
    require('actual_sscard_url_in_manifest', sscard['url'] == sscard['resolved_url'] == expected_url)
    require('actual_sscard_bytes_hash', sscard['bytes'] == 2279120 and sscard['sha256']
            == '63c8342280dee3bd67c1322141ef42bc6a350e479cbaeb69fc4474a67b14796d')
    require('assay_metadata_unit_semantics', spec['sscard']['unit'] == {'SSB2M':'mg/L','SSCYST':'mg/L'}
            and spec['sscard']['matrix'] == 'serum'
            and spec['sscard']['platform'] == 'Siemens Dimension Vista 1500'
            and spec['sscard']['last_revised'] == '2022-07'
            and spec['sscard']['detection_limits_mg_l'] == {'SSB2M':[0.72,23.0],'SSCYST':[0.23,8.0]}
            and spec['sscard']['beta_trace_protein_present'] is False)
    require('two_exact_pdf_denials_documented', manifest['documentation_not_retrieved']
            == a['documentation_not_retrieved'] and len(manifest['documentation_not_retrieved']) == 2)
    for row in manifest['documentation_not_retrieved']:
        require('not_retried_pdf:' + row['file'], row['status'] == 'denied_not_retried'
                and row['original_http_status'] == 403 and row['original_run_id'] == 40907895237
                and row['original_job_id'] == 123567576780 and row['file'] not in rows
                and 'sha256' not in row and 'bytes' not in row)
    require('exact_denied_pdf_url_set', {r['url'] for r in manifest['documentation_not_retrieved']} == {
        spec['mortality_dictionary_url'], spec['mortality_description_url']})

    tree = ast.parse((folder/'cohort.py').read_text())
    fn = next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='six_year_weights')
    actual = [n for n in ast.walk(fn) if isinstance(n,ast.Call)
              and isinstance(n.func,ast.Attribute) and isinstance(n.func.value,ast.Name)
              and n.func.value.id=='np' and n.func.attr=='where']
    expected = ast.parse('np.where(frame.SDDSRVYR.isin([1,2]),(2/3)*frame.WTSSCB4Y,(1/3)*frame.WTSSCB2Y)',mode='eval').body
    require('weight_formula_static_ast_exact', len(actual)==1
            and ast.dump(actual[0],include_attributes=False)==ast.dump(expected,include_attributes=False))
    require('full_sscard_row_count', a['corrected_sscard_source_rows']==21890)
    zero = {r['filename']:r for r in a['xpt_zero_audits']}
    require('twentytwo_xpt_zero_audits', len(zero)==len(a['xpt_zero_audits'])==22
            and set(zero)=={n for n in rows if n.endswith('.xpt')})
    for name, z in zero.items():
        require('xpt_zero_provenance:' + name, z['sha256']==rows[name]['sha256']
                and z['epsilon_threshold_used'] is False)
    require('exported_sscard_flag0_and_pristine0_counts', zero['1999_SSCARD_A.xpt']['nrows']==21890
            and zero['1999_SSCARD_A.xpt']['numeric_zero_counts']['SSB2ML']==21824
            and zero['1999_SSCARD_A.xpt']['numeric_zero_counts']['SSCYSTL']==21824
            and zero['1999_SSCARD_A.xpt']['numeric_zero_counts']['SSPRIS']==1674)
    readers = a['mortality_reader_audits']
    require('three_mortality_readers_exact_export', len(readers)==3
            and {r['cycle'] for r in readers}=={1999,2001,2003})
    for r in readers:
        require('mortality_reader:' + str(r['cycle']), r['exact_match'] is True and r['fields']==8
                and r['cell_comparisons']==r['rows']*8 and r['sha256']==rows[r['file']]['sha256'])
    require('mortality_total_rows_and_cells', sum(r['rows'] for r in readers)==31126
            and sum(r['cell_comparisons'] for r in readers)==249008)

    full, eligible = a['full_DEMO_masked_design'], a['eligible_masked_design']
    require('exact_full_eligible_inventory_equality', full==eligible and a['full_DEMO_design_inventory_preserved'] is True)
    inv = {int(s):list(p) for s,p in full.items()}
    require('fortythree_exact_masked_strata', set(inv)==set(range(1,44)))
    require('masked_psu_codes_valid', all(ps==sorted(set(ps)) and len(ps)>=2
            and set(ps)<=({1,2,3} if st<=13 else {1,2}) for st,ps in inv.items()))
    distribution = dict(Counter(len(ps) for ps in inv.values()))
    psu_n = sum(len(ps) for ps in inv.values())
    require('eightyseven_psu_distribution_including_m3', psu_n==87 and distribution=={3:1,2:42}
            and a['eligible_psu_count_distribution']=={'2':42,'3':1}
            and a['full_DEMO_masked_psu_n']==a['eligible_masked_psu_n']==87 and inv[1]==[1,2,3])

    flow = a['cohort_flow']
    require('cycle_demographic_counts', [(r['cycle'],r['survey_n']) for r in flow]
            == [(1999,9965),(2001,11039),(2003,10122)])
    for r in flow:
        values=[r[k] for k in ('survey_n','age40_75_n','MEC_target_n','linked_target_n','eligible_weighted_subsample_n')]
        require('nested_flow:' + str(r['cycle']), all(isinstance(x,int) and x>=0 for x in values)
                and values==sorted(values,reverse=True)
                and r['linked_target_n']-r['eligible_weighted_subsample_n']==r['excluded_no_corrected_subsample_weight_n'])
    require('eligible_flow_sum', sum(r['eligible_weighted_subsample_n'] for r in flow)==a['eligible_n']==results['eligible_n']==6546)
    primary,sens=a['common_domain'],a['reported_sensitivity_domain']
    require('primary_target_counts', (primary['n'],primary['events8'])==(6044,586)
            and results['n']==6044 and results['events8']==586 and results['strata_n']==43 and results['psu_n']==87)
    require('sensitivity_target_counts', (sens['n'],sens['events8'])==(6054,586)
            and results['reported_values_sensitivity']['n']==6054 and results['reported_values_sensitivity']['events8']==586)
    for label,dom,flags in [('primary',primary,{'0.0'}),('sensitivity',sens,{'0.0','1.0'})]:
        require('domain_sums:' + label, sum(dom['cycles'].values())==dom['n']
                and sum(dom['pristine_counts'].values())==dom['n'] and dom['events8']<=dom['n']<=6546
                and dom['weight_sum']>0 and 0<dom['effective_n']<=dom['n'])
        for flag,v in dom['flag_counts'].items():
            require('domain_flag_counts:' + label+':'+flag, set(v)<=flags and sum(v.values())==dom['n'])
    require('coverage_root_audit_equal', 0<a['common_weighted_coverage']<=1
            and a['common_weighted_coverage']==results['common_weighted_coverage'])
    require('effective_n_root_audit_equal', close(primary['effective_n'],results['effective_n']))
    patterns=a['exclusion_reason_patterns']; reason_names=set(a['exclusion_reason_counts'])
    require('exclusion_pattern_total', all(isinstance(n,int) and n>=0 for n in patterns.values())
            and sum(patterns.values())==6546 and patterns['complete']==6044)
    reconstructed={r:0 for r in reason_names};multiple=0
    for label,n in patterns.items():
        if label=='complete':continue
        reasons=label.split('|')
        require('known_unique_pattern_reasons:' + label, len(reasons)==len(set(reasons)) and set(reasons)<=reason_names)
        for r in reasons:reconstructed[r]+=n
        if len(reasons)>1:multiple+=n
    require('marginal_exclusion_reasons_reconstructed', reconstructed==a['exclusion_reason_counts'])
    require('survivor_followup_gate', a['minimum_survivor_followup_months']==177 and a['minimum_survivor_followup_months']>=96)
    require('zero_month_events_disclosed', a['zero_month_deaths']==0)
    models={'routine','routine_b2m','routine_cys','routine_cys_b2m'}
    require('four_primary_arms_one_domain_counts', set(results['point'])==models and all(
        p['n']==6044 and p['events']==586 and p['n_positive_weight']==6044 and p['events_positive_weight']==586
        and close(p['weight_sum'],primary['weight_sum']) for p in results['point'].values()))
    require('four_sensitivity_arms_one_domain_counts', set(results['reported_values_sensitivity']['point'])==models and all(
        p['n']==6054 and p['events']==586 and p['n_positive_weight']==6054 and p['events_positive_weight']==586
        and close(p['weight_sum'],sens['weight_sum']) for p in results['reported_values_sensitivity']['point'].values()))
    require('same_common_domain_assertion', results['same_common_domain_all_four_models'] is True)
    require('forty_fold_records', len(folds)==40)
    for domain, n, events in [('primary_common',6044,586),('reported_values_sensitivity',6054,586)]:
        group=[r for r in folds if r['analysis_domain']==domain]
        # Use actual exported domain labels without assuming identity from per-arm counts alone.
        if not group and domain=='reported_values_sensitivity':
            other={r['analysis_domain'] for r in folds}-{'primary_common'}
            require('one_sensitivity_fold_domain',len(other)==1)
            group=[r for r in folds if r['analysis_domain'] in other]
        require('twenty_domain_fold_records:'+domain,len(group)==20)
        for model in models:
            entries=sorted([r for r in group if r['model']==model],key=lambda r:r['outer_fold'])
            require('five_scoring_folds:'+domain+model,len(entries)==5 and [r['outer_fold'] for r in entries]==list(range(5))
                    and sum(r['scoring_n'] for r in entries)==n and sum(r['scoring_events'] for r in entries)==events
                    and sum(r['scoring_psu_n'] for r in entries)==87)
            for e in entries:
                require('fold_domain_arithmetic:'+domain+model+str(e['outer_fold']),
                        e['training_n']+e['scoring_n']==n and e['training_events']+e['scoring_events']==events
                        and e['training_psu_n']+e['scoring_psu_n']==87 and e['PSU_overlap_n']==0)
        for fold in range(5):
            entries=[r for r in group if r['outer_fold']==fold]
            keys=['training_n','scoring_n','training_events','scoring_events','training_psu_n','scoring_psu_n','training_weight_sum']
            require('matching_arm_fold_aggregates:'+domain+str(fold),len({tuple(r[k] for k in keys) for r in entries})==1)
    require('no_independence_or_clinical_claims', results['independent_external_validation'] is False
            and results['clinical_use_ready'] is False and results['participant_records_exported'] is False
            and results['previous_frozen_coefficients_used'] is False)
    require('reuse_history_disclosed', a['source_all_participants_from_previously_loaded_release'] is True
            and a['prior_Cys_fitted'] is True and a['prior_B2M_loaded'] is True
            and a['old_fitted_M4_participant_intersection_not_reconstructed'] is True)

    require('log_checkout_merge_and_science_parent', ('HEAD is now at a752d10 Merge '+CODE+' into ') in log and MERGE in log)
    require('log_fifty_source_success', '{"retrieved": 50, "failures": []}' in log)
    point_line=next(l for l in log.splitlines() if ' B2M31_POINT_COMPUTED ' in l)
    point_time=dt.datetime.fromisoformat(point_line.split(' ',1)[0].replace('Z','+00:00'))
    commit_time=dt.datetime.fromisoformat(subprocess.check_output(['git','show','-s','--format=%cI',CODE],cwd=repo,text=True).strip())
    require('amendment_before_comparison',commit_time<point_time)
    require('amendment_pinned', 'SOURCE_INTAKE_AMENDMENT.md' in metadata)
    def marked(begin,end):
        lines=[re.sub(r'^\d{4}-\d\d-\d\dT[^ ]+ ', '',line) for line in log.splitlines()]
        i=lines.index(begin);j=lines.index(end)
        return strict('\n'.join(lines[i+1:j]))
    require('result_equals_marked_job_export',marked('B2M31_RESULTS_BEGIN','B2M31_RESULTS_END')==results)
    require('manifest_equals_marked_job_export',marked('B2M31_MANIFEST_BEGIN','B2M31_MANIFEST_END')==manifest)
    require('folds_equal_marked_job_export',marked('B2M31_PARAMETERS_BEGIN','B2M31_PARAMETERS_END')==folds)
    require('calibration_equal_marked_job_export',marked('B2M31_CALIBRATION_BINS_BEGIN','B2M31_CALIBRATION_BINS_END')==calibration)
    for name,key in [('bootstrap_aggregates.csv','bootstrap_aggregates_csv_sha256'),
                     ('primary_bootstrap.csv','primary_bootstrap_csv_sha256'),
                     ('secondary_bootstrap.csv','secondary_bootstrap_csv_sha256')]:
        require('aggregate_csv_bytes_hash:'+name,digest((folder/name).read_bytes())==results[key])

    return {
        'schema':'b2m31-independent-source-audit/1','status':'passed_with_explicit_scope_limits',
        'audit_date':'2026-10-29','no_source_records_opened':True,'no_network_requests':True,'no_model_refits':True,
        'scientific_head':CODE,'checked_out_PR_merge_commit_from_log':MERGE,
        'run_id_supplied_by_parent':40908831254,'job_id_supplied_by_parent':123570909752,
        'script_sha256':digest(Path(__file__).read_bytes()),
        'artifacts':{p.name:{'sha256':digest(p.read_bytes()),'bytes':p.stat().st_size} for p in
                     [folder/'RESULTS_SUMMARY.json',folder/'SOURCE_MANIFEST.json',folder/'FOLD_PARAMETERS.json',folder/'CALIBRATION_BINS.json',log_path]},
        'checks':checks,'passed_check_count':len(checks),'failed_check_count':0,
        'locked_metadata_count':len(metadata),'locked_historical_python_count':len(historical_code),
        'scientific_export_hashes_checked_against_worktree_and_scientific_head':len(scientific),
        'manifest':{'files':50,'historical_data_files':25,'full_source_bytes_reported':sum(r['bytes'] for r in rows.values()),
                    'data_source_bytes_reported':sum(r['bytes'] for r in data.values()),
                    'manifest_sha256':digest(manifest_bytes),'source_content_sha256':content_hash,
                    'actual_sscard_entry':sscard,'documentation_denials':manifest['documentation_not_retrieved']},
        'assay':{'matrix':'serum','reported_units':'mg/L','platform':'Siemens Dimension Vista 1500','release':'2022-07',
                 'full_marginal_policy':'code guards enforced; exported flag0/SSPRIS0 zero-count audit independently crosschecked; no raw counts recomputed',
                 'LOD_not_LOQ':True,'BTP_present':False},
        'weight_formula_ast_verified':{'SDDSRVYR_1_2':'(2/3)*WTSSCB4Y','SDDSRVYR_3':'(1/3)*WTSSCB2Y'},
        'design':{'strata_n':43,'PSU_n':87,'PSU_per_stratum_distribution':distribution,'m3_strata':{str(k):v for k,v in inv.items() if len(v)==3},
                  'full_eligible_inventory_exactly_equal':True,'primary_common_has_87_clusters_from_fold_aggregates':True},
        'domain':{'eligible_n':6546,'primary':primary,'sensitivity':sens,'flow':flow,'excluded_union_n':6546-6044,
                  'marginal_exclusion_sum':sum(reconstructed.values()),'multiple_exclusion_reason_people_n':multiple,
                  'sensitivity_added_n':10,'sensitivity_added_events':0,
                  'sensitivity_added_by_cycle':{c:sens['cycles'][c]-primary['cycles'][c] for c in primary['cycles']},
                  'sensitivity_added_pristine':{k:sens['pristine_counts'].get(k,0)-primary['pristine_counts'].get(k,0) for k in sens['pristine_counts']},
                  'minimum_survivor_months':a['minimum_survivor_followup_months']},
        'chronology':{'amendment_commit_time':commit_time.isoformat(),'first_comparative_point_log_time':point_time.isoformat(),
                      'amendment_precedes_comparison':True},
        'limitations':[
            'No participant/raw-source bytes were opened or independently rehashed; fifty fingerprints were checked as manifest evidence and twentyfive against exact historical pins.',
            'Runtime SSCARD anchor resolution is bound by pinned resolver and completed source gate; raw index HTML was not independently re-parsed in this audit.',
            'Full flag1/flag2/missing/pristine1 marginal categories are enforced by pinned validation code, not all explicitly exported for independent numerical recomputation.',
            'Cycle-specific unscaled/scaled weight sums and eligible total weight are not exported; per-row W6 and weighted-coverage denominator were not independently recomputed.',
            'Counts and shared-fold aggregates cannot independently prove participant-level mask equality; pinned code applies one common mask/fold array to all four arms.',
            'Custom refit-bootstrap numerical intervals are audited separately; masked design preservation alone does not establish official survey-design coverage.',
            'Reused SSCARD_A participants and previous Cys fitting/B2M loading; exact intersection with former fitted M4 domain not reconstructed; not independent external validation.'
        ]}


def markdown(doc):
    m = doc['manifest']
    d = doc['domain']
    a = m['actual_sscard_entry']
    return f"""# B2M31: независимый аудит источников и домена

**Пройден с явно ограниченной областью.** Выполнено {doc['passed_check_count']} проверок JSON, хешей, кода и агрегатной арифметики; отказов — 0. Исходные записи и байты наборов не открывались. Сетевых запросов, загрузки ZIP и повторного обучения не было. Численные интервалы bootstrap проверяет отдельный reviewer.

В журнале checkout указан PR merge `{MERGE}`, содержащий scientific head `{CODE}`. Все 90 экспортированных хешей научных файлов совпали с этим head и текущими файлами. Проверены 12 хешей метаданных LOCK и {doc['locked_historical_python_count']} хеша прежнего Python-кода. Amendment закреплён {doc['chronology']['amendment_commit_time']}; первая сравнительная точка в журнале — {doc['chronology']['first_comparative_point_log_time']}. Поправка предшествует сравнению.

Множество 50 источников точно совпало с независимо сохранённым reference. У всех записей — retrieved, официальный HTTPS URL, положительный ограниченный byte-length и SHA256. Сумма заявленных байтов — {m['full_source_bytes_reported']:,}; для 25 файлов данных — {m['data_source_bytes_reported']:,}. Все 25 data SHA точно совпали по имени с прежними pins nhanes04. SHA исходных байтов экспортированного SOURCE_MANIFEST: `{m['manifest_sha256']}`; canonical content SHA: `{m['source_content_sha256']}`. Это проверка манифеста, а не независимое повторное хеширование самих наборов.

Фактический SSCARD URL из runtime manifest: `{a['url']}`. Resolved URL совпадает; размер — {a['bytes']:,} байтов; SHA — `{a['sha256']}`. Закреплённый resolver требует уникальный официальный anchor и запрещает угадывание URL; исходный index HTML я повторно не разбирал. Serum, mg/L, Dimension Vista 1500 и релиз July 2022 согласованы с SOURCE_SPEC. Число 21890 и экспортированные счётчики нулей двух flag0 — по 21824, SSPRIS0 — 1674 — проверены. Полные категории flag1/flag2/missing проходят закреплённые кодовые проверки, но не все отдельно экспортированы. Независимый подсчёт этих категорий по исходным записям не заявляется. LOD не назван LOQ; BTP отсутствует.

Два mortality PDF с прежним HTTP403 сохранены отдельно как `denied_not_retried`, с run 40907895237 / job 123567576780. Они не входят в 50 успешных источников и не имеют фиктивных новых SHA/байтов. Три mortality-reader audit дают 31126 строк / 249008 сопоставлений ячеек с точным совпадением. Новое скачивание этих PDF не заявлено.

Формула W6 проверена через AST: 2/3 × WTSSCB4Y для 1999–2002 и 1/3 × WTSSCB2Y для 2003–2004. Поцикловые суммы исходных и масштабированных весов, а также полный знаменатель eligible weight не экспортированы. Поэтому применение множителей к каждой записи и знаменатель coverage независимо не пересчитаны.

Полный и eligible masked design совпадают по содержимому: 43 strata / 87 PSU. В 42 strata по две PSU; в stratum 1 — PSU [1, 2, 3]. Две PSU для всех strata не предполагаются. Агрегаты фолдов подтверждают 87 common clusters, нулевой заявленный PSU overlap и одинаковые счётчики всех четырёх моделей.

| Домен | n | Смертей за 8 лет | 1999 / 2001 / 2003 | Non-pristine / pristine |
|---|---:|---:|---|---|
| Eligible | 6546 | Не экспортировано совокупно | 1993 / 2330 / 2223 | Не экспортировано |
| Primary, оба flags0 | 6044 | 586 | 1889 / 2131 / 2024 | 43 / 6001 |
| Reported, flags0/1 | 6054 | 586 | 1891 / 2132 / 2031 | 43 / 6011 |

Исключено 502 человека. Сумма отдельных причин — 530 из-за перекрытия; у 24 человек несколько причин. Все marginal reason counts независимо восстановлены из patterns; сумма patterns, включая complete 6044, равна eligible 6546. В primary оба flags равны 0; sensitivity содержит только 0/1. Sensitivity добавляет 10 человек: по циклам +2 / +1 / +7, все pristine, дополнительных смертей за 8 лет — 0. Абсолютные scores между разными доменами не объявляются улучшением. Минимальный survivor follow-up — 177 месяцев, что превышает 96; смертей на нулевом месяце — 0.

Сумма весов primary — {d['primary']['weight_sum']:.8f}, effective n — {d['primary']['effective_n']:.8f}, coverage — 0.9294070995654213. Проверена согласованность RESULT, point, source audit и всех четырёх моделей. Отсутствующий знаменатель не подменён реконструкцией из округлённого процента.

SSCARD_A и участники использовались ранее: Cys уже оценён в M4, B2M загружался. Точное пересечение с прежним fitted M4 domain не восстановлено. Итерация 31 — exploratory B2M contrast после III28/29, не новая когорта, untouched holdout или независимая валидация freeze30. Проверены закреплённые хеши FROZEN_MODELS/FREEZE_MANIFEST и прежнего научного кода. Пределы независимой проверки перечислены в companion JSON.
"""


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--repo',type=Path,default=Path('/workspace/aging_theories'))
    p.add_argument('--log',type=Path,default=Path('/tmp/b2m31_success_job.log'))
    p.add_argument('--reference',type=Path,default=Path('/tmp/b2m31_expected_manifest_reference.json'))
    p.add_argument('--out',type=Path,default=Path('/tmp/b2m31_source_audit'))
    args=p.parse_args();doc=audit(args.repo,args.log,args.reference)
    args.out.with_suffix('.json').write_text(json.dumps(doc,ensure_ascii=False,indent=2,allow_nan=False)+'\n')
    args.out.with_suffix('.md').write_text(markdown(doc))
    print(json.dumps({'status':doc['status'],'checks':doc['passed_check_count'],'failures':0,
                      'files':50,'historical_data_files':25,'strata':43,'PSU':87,'n':6044,'events':586}))
