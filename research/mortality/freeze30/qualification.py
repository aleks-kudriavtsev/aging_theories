"""Pure metadata qualification of a candidate external frozen-model evaluation.

No file/network/participant operations occur here. Declarations and URLs do not
prove legal permission, independence, assay interchangeability or chronology.
The reviewed external-mapping registry is intentionally EMPTY: metadata alone
cannot authorize empirical scoring, including a fully populated example.

Schema ``external_qualification30/1`` (see ``blank_contract``):
- cohort_id/name/release_id, public_unrestricted data_class, execution_context;
- sources: official_url, publisher_domain, release_id, kind=person_level,
  direct_download=True, requires_account=False, requires_agreement=False;
- overlap: checked_against list of known programme cohorts, disjoint_evidence_urls,
  assessment=independent_participant_universe_documented;
- endpoint: baseline_exam origin, all_cause_death, target ages40..75, horizon8,
  explicit time unit/conversion, linkage, verified survivor followup and censors;
- routine_core:11 engineered terms with exact semantic_id/unit, source variables,
  raw_inputs, definition_reference, baseline_exam and derivation_reviewed=True;
- raw_measurements: exact units/matrices, source/visit/method and missing policy;
- proteins:CST3/B2M/BTP with exact identity, mg/L serum absolute concentration,
  method/version/quantification policy and documented historical transport;
- sampling: documented complex probability or explicitly unweighted cohort;
- protocol:SHA256, target bundleSHA256, immutable lock reference, pre-outcome dates;
- mapping_id: a future code-reviewed registry key, never self-approval.

Historical NHANESIII CEP calibration cannot be blindly applied to modern IDMS
creatinine. Routine/raw and protein transport claims need reviewed mapping code.
No unknown, controlled/account or short-followup path is approved by this module.
"""
from __future__ import annotations
import datetime as _dt
import re
import unicodedata
from urllib.parse import urlsplit

SCHEMA = 'external_qualification30/1'
QUALIFIED_EXTERNAL_MAPPINGS = {}  # Requires a separately reviewed code commit.
REUSED_PROGRAMME_COHORTS = ('NHANESIII', 'NHANES1999_2018', 'NHEFS', 'FLCHAIN', 'PBC')
CORE_FEATURES = {
    'age_linear': ('interview_age_minus60_div10', 'dimensionless', ('age_years',)),
    'age_quadratic': ('square_of_age_linear', 'dimensionless', ('age_years',)),
    'male': ('recorded_male1_female0', 'dimensionless', ('sex',)),
    'current': ('current_cigarettes_after100_lifetime', 'dimensionless', ('ever100_cigarettes','current_cigarettes')),
    'former': ('former_cigarettes_after100_lifetime', 'dimensionless', ('ever100_cigarettes','current_cigarettes')),
    'sbp': ('baseline_mean_sbp_minus130_div20', 'dimensionless', ('systolic_bp',)),
    'bmi': ('baseline_bmi_minus27_div5', 'dimensionless', ('bmi',)),
    'hba1c': ('baseline_DCCT_aligned_hba1c_percent', '%', ('hba1c',)),
    'egfr_low': ('CKDEPI2021_creatinine_min60_div15', 'dimensionless', ('serum_creatinine','age_years','sex')),
    'egfr_high': ('CKDEPI2021_creatinine_above60_div30', 'dimensionless', ('serum_creatinine','age_years','sex')),
    'log_uacr': ('log2_UACR_mg_per_g', 'dimensionless', ('urine_albumin','urine_creatinine')),
}
RAW_INPUTS = {
    'age_years': ('years','not_applicable'), 'sex': ('coded','not_applicable'),
    'ever100_cigarettes': ('coded','not_applicable'), 'current_cigarettes': ('coded','not_applicable'),
    'systolic_bp': ('mmHg','not_applicable'), 'bmi': ('kg/m2','not_applicable'),
    'hba1c': ('%','whole_blood'), 'serum_creatinine': ('mg/dL','serum'),
    'urine_albumin': ('ug/mL','urine'), 'urine_creatinine': ('mg/dL','urine'),
}
PROTEINS = {'CST3':'P01034', 'B2M':'P61769', 'BTP':'P41222'}
HISTORICAL_ASSAYS = {
    'CST3':'NHANESIII_original2006_DadeBehring_NLatex_PENIA',
    'B2M':'NHANESIII_original2009_Siemens_NLatex_B2M',
    'BTP':'NHANESIII_original2009_Siemens_NLatex_BTP',
}
_TOP = {'schema_version','cohort_id','cohort_name','release_id','data_class','execution_context',
        'sources','overlap','endpoint','routine_core','raw_measurements','proteins','sampling','protocol','mapping_id'}
_CONTEXTS = {'ai_session','github_actions','authorized_human_local'}
_CATEGORIES = ('access','overlap','endpoint','routine','assay','sampling','protocol','mapping')


def _text(value):
    return (isinstance(value,str) and bool(value.strip())
            and value.strip().lower() not in {'unknown','none','pending','not_reported','unresolved'})


def _url(value):
    if not _text(value): return False
    try:
        p=urlsplit(value)
        return (p.scheme=='https' and bool(p.hostname) and p.username is None and p.password is None
                and p.port in (None,443) and bool(p.path) and not p.fragment)
    except (ValueError,TypeError): return False


def _sha(value):
    return isinstance(value,str) and re.fullmatch(r'[0-9a-f]{64}',value) is not None


def _normal(value):
    return re.sub('[^A-Z0-9]','',unicodedata.normalize('NFKC',value).upper()) if isinstance(value,str) else ''


def _obj(c,key):
    value=c.get(key,{})
    if not isinstance(value,dict): raise ValueError(key+'_must_be_metadata_object')
    return value


def _instant(value):
    if not isinstance(value,str): return None
    try:
        result=_dt.datetime.fromisoformat(value.replace('Z','+00:00'))
        return result if result.tzinfo is not None else None
    except ValueError: return None


def blank_contract():
    """Return an unqualified metadata template; no participant fields or IO."""
    routine={name:dict(available=False,semantic_id=semantic,unit=unit,raw_inputs=list(raw),
                       source_variables=[],definition_reference=None,baseline_visit=None,derivation_reviewed=False)
             for name,(semantic,unit,raw) in CORE_FEATURES.items()}
    raw={name:dict(available=False,source_variable=None,unit=unit,matrix=matrix,baseline_visit=None,
                   method=None,method_version=None,definition_reference=None,missing_and_substitution_policy=None,
                   calibration_basis=None,transport={}) for name,(unit,matrix) in RAW_INPUTS.items()}
    proteins={name:dict(available=False,source_variable=None,uniprot=identity,unit='mg/L',matrix='serum',
                        scale='absolute_concentration',baseline_visit=None,method=None,method_version=None,
                        quantification_policy=None,source_reference=None,transport={})
              for name,identity in PROTEINS.items()}
    return dict(schema_version=SCHEMA,cohort_id=None,cohort_name=None,release_id=None,
        data_class=None,execution_context=None,sources=[],overlap={},endpoint={},routine_core=routine,
        raw_measurements=raw,proteins=proteins,sampling={},protocol={},mapping_id=None)


def _transport(row, expected, prefix, problems):
    t=row.get('transport',{})
    if not isinstance(t,dict): raise ValueError(prefix+'_transport_must_be_metadata_object')
    kind=t.get('kind')
    if kind=='documented_same_historical_method':
        if t.get('historical_method_id')!=expected or not _url(t.get('evidence_url')):
            problems.append(prefix+':historical_method_equivalence_unresolved')
    elif kind=='paired_bridge_to_historical_method':
        if (t.get('historical_method_id')!=expected or not _url(t.get('evidence_url'))
                or not _text(t.get('bridge_id'))):
            problems.append(prefix+':historical_bridge_specification_unresolved')
        # A claimed paired bridge is not a validated implementation.
        problems.append(prefix+':paired_bridge_not_in_reviewed_external_mapping')
    else:
        problems.append(prefix+':historical_assay_transport_unresolved')


def audit(c):
    """Return blockers only; never authenticate rights or enable unreviewed IO."""
    if not isinstance(c,dict): raise ValueError('metadata_object_required')
    if set(c)-_TOP: raise ValueError('unsupported_top_level_metadata_fields')
    b={key:[] for key in _CATEGORIES}
    if c.get('schema_version')!=SCHEMA: b['access'].append('unsupported_schema')
    if not _text(c.get('cohort_id')) or not _text(c.get('release_id')): b['access'].append('exact_cohort_and_release_required')
    if c.get('data_class')!='public_unrestricted': b['access'].append('only_direct_public_unrestricted_data_supported')
    context=c.get('execution_context')
    if not isinstance(context,str) or context not in _CONTEXTS: b['access'].append('execution_context_unknown_or_unsupported')
    identity=_normal(c.get('cohort_id'))+_normal(c.get('cohort_name'))
    sources=c.get('sources',[])
    if not isinstance(sources,list): raise ValueError('sources_must_be_metadata_list')
    if not sources: b['access'].append('direct_official_person_files_required')
    hrs=('HRS' in identity or 'HEALTHANDRETIREMENTSTUDY' in identity)
    for i,row in enumerate(sources):
        if not isinstance(row,dict): raise ValueError('source_must_be_metadata_object')
        url=row.get('official_url');domain=row.get('publisher_domain')
        host=urlsplit(url).hostname if _url(url) else None
        if host and (host=='umich.edu' or host.endswith('.umich.edu')): hrs=True
        if (not host or not _text(domain) or not (host==domain or host.endswith('.'+domain))
                or row.get('release_id')!=c.get('release_id') or row.get('kind')!='person_level'):
            b['access'].append(f'source{i}:official_release_person_file_unresolved')
        if row.get('direct_download') is not True: b['access'].append(f'source{i}:direct_download_unconfirmed')
        for flag in ('requires_account','requires_agreement'):
            if row.get(flag) is not False: b['access'].append(f'source{i}:{flag}_or_unknown')
        if not _url(row.get('public_release_reference')): b['access'].append(f'source{i}:public_release_evidence_required')
    if hrs and context in {'ai_session','github_actions'}:
        b['access'].append('HRS_person_records_prohibited_in_AI_or_unattended_cloud')

    overlap=_obj(c,'overlap')
    if any(x in identity for x in ('NHANESIII','NHANES3')):
        b['overlap'].append('NHANESIII_first_second_or_renamed_release_already_used')
    if 'NHANES' in identity and any(str(year) in identity for year in range(1999,2019)):
        b['overlap'].append('NHANES1999_2018_programme_cycle_already_used')
    if any(identity.startswith(name) for name in ('NHEFS','FLCHAIN','PBC')):
        b['overlap'].append('known_previously_used_programme_cohort')
    checked=overlap.get('checked_against')
    refs=overlap.get('disjoint_evidence_urls')
    if (overlap.get('assessment')!='independent_participant_universe_documented'
            or not isinstance(checked,list) or not all(isinstance(x,str) for x in checked)
            or not set(REUSED_PROGRAMME_COHORTS)<=set(checked)
            or not isinstance(refs,list) or not refs or not all(_url(u) for u in refs)):
        b['overlap'].append('documented_disjoint_participant_universe_not_established')

    e=_obj(c,'endpoint')
    for field,expected in [('kind','all_cause_death'),('baseline_origin','baseline_exam'),
                           ('target_age_min_years',40),('target_age_max_years',75),('horizon_years',8),
                           ('baseline_visit','baseline_exam')]:
        actual=e.get(field)
        if actual!=expected or (isinstance(expected,int) and type(actual) not in (int,float)):
            b['endpoint'].append(field+'_not_canonical')
    if e.get('linkage_same_person_verified') is not True: b['endpoint'].append('same_person_linkage_unverified')
    if not _url(e.get('definition_reference')): b['endpoint'].append('endpoint_definition_evidence_required')
    unit=e.get('followup_unit');conversion=e.get('conversion_to_years')
    if (not isinstance(unit,str) or not isinstance(conversion,str)
            or (unit,conversion) not in {('years','identity'),('months','divide_by_12')}):
        b['endpoint'].append('followup_unit_or_conversion_unresolved')
    minimum=e.get('minimum_survivor_followup_years')
    if type(minimum) not in (int,float) or not 8<=minimum<1000:
        b['endpoint'].append('eight_year_survivor_followup_not_established')
    if type(e.get('early_censored_n')) is not int or e['early_censored_n']!=0:
        b['endpoint'].append('early_censoring_requires_separate_IPCW_protocol_and_implementation')
    if e.get('followup_support_verified') is not True: b['endpoint'].append('followup_support_not_verified')

    raw=_obj(c,'raw_measurements');core=_obj(c,'routine_core')
    for name,(semantic,unit,inputs) in CORE_FEATURES.items():
        row=core.get(name)
        if not isinstance(row,dict): b['routine'].append(name+':missing_engineered_term');continue
        if row.get('available') is not True or row.get('semantic_id')!=semantic or row.get('unit')!=unit:
            b['routine'].append(name+':definition_unit_or_availability_unresolved')
        declared=row.get('raw_inputs');variables=row.get('source_variables')
        if (not isinstance(declared,list) or not all(isinstance(x,str) for x in declared) or set(declared)!=set(inputs)
                or not isinstance(variables,list) or not variables or not all(_text(x) for x in variables)):
            b['routine'].append(name+':exact_raw_mapping_unresolved')
        expected_sources=[raw.get(key,{}).get('source_variable') if isinstance(raw.get(key),dict) else None for key in inputs]
        if (not all(_text(x) for x in expected_sources) or not isinstance(variables,list)
                or not all(_text(x) for x in variables) or set(variables)!=set(expected_sources)):
            b['routine'].append(name+':engineered_to_raw_source_variables_disagree')
        if row.get('baseline_visit')!='baseline_exam' or row.get('derivation_reviewed') is not True or not _url(row.get('definition_reference')):
            b['routine'].append(name+':baseline_or_derivation_evidence_unresolved')
    for name,(unit,matrix) in RAW_INPUTS.items():
        row=raw.get(name)
        if not isinstance(row,dict): b['routine'].append(name+':raw_source_missing');continue
        if (row.get('available') is not True or row.get('unit')!=unit or row.get('matrix')!=matrix
                or row.get('baseline_visit')!='baseline_exam' or not _text(row.get('source_variable'))
                or not _url(row.get('definition_reference')) or not _text(row.get('missing_and_substitution_policy'))):
            b['routine'].append(name+':raw_definition_unit_visit_or_missing_policy_unresolved')
        if name in ('hba1c','serum_creatinine','urine_albumin','urine_creatinine'):
            if not _text(row.get('method')) or not _text(row.get('method_version')):
                b['assay'].append(name+':method_version_unresolved')
        if name=='serum_creatinine':
            if row.get('calibration_basis')!='NHANESIII_original_CEP':
                b['assay'].append('serum_creatinine:modern_IDMS_or_other_scale_cannot_use_historical_CEP_equation_without_reviewed_bridge')
            _transport(row,'NHANESIII_original_CEP_before_0.960x_minus0.184',name,b['assay'])
    proteins=_obj(c,'proteins')
    for name,uniprot in PROTEINS.items():
        row=proteins.get(name)
        if not isinstance(row,dict): b['assay'].append(name+':protein_missing');continue
        for field,expected in [('uniprot',uniprot),('unit','mg/L'),('matrix','serum'),('scale','absolute_concentration'),('baseline_visit','baseline_exam')]:
            if row.get(field)!=expected: b['assay'].append(name+':'+field+'_mismatch')
        for field in ('source_variable','method','method_version','quantification_policy'):
            if not _text(row.get(field)): b['assay'].append(name+':'+field+'_unresolved')
        if row.get('available') is not True or not _url(row.get('source_reference')):
            b['assay'].append(name+':availability_or_primary_evidence_unresolved')
        _transport(row,HISTORICAL_ASSAYS[name],name,b['assay'])

    sampling=_obj(c,'sampling')
    if sampling.get('documented') is not True or not _url(sampling.get('design_reference')):
        b['sampling'].append('sampling_and_target_population_not_documented')
    if sampling.get('same_joint_domain_supported') is not True: b['sampling'].append('joint_panel_weight_support_unverified')
    design=sampling.get('design')
    if design=='complex_probability':
        if not all(_text(sampling.get(k)) for k in ('weight_variable','stratum_variable','PSU_variable')):
            b['sampling'].append('weights_stratum_PSU_mapping_required')
    elif design=='unweighted_cohort':
        if sampling.get('analysis_weight_rule')!='unit_weights_received_cohort_estimand_only':
            b['sampling'].append('unit_weights_require_explicit_received_cohort_estimand')
    else: b['sampling'].append('sampling_design_unresolved')
    if sampling.get('weight_variable')=='WTCYPEX6' and 'NHANES' not in identity:
        b['sampling'].append('NHANES_specific_subsample_weight_cannot_be_external_default')

    p=_obj(c,'protocol')
    if not _sha(p.get('sha256')) or not _sha(p.get('target_bundle_sha256')):
        b['protocol'].append('protocol_and_exact_target_bundle_sha256_required')
    if not _url(p.get('immutable_lock_reference')): b['protocol'].append('immutable_protocol_lock_reference_required')
    frozen=_instant(p.get('frozen_at_utc'));outcome=_instant(p.get('first_outcome_access_at_utc'))
    accessed=p.get('outcomes_accessed')
    if p.get('fixed_before_outcomes') is not True or frozen is None:
        b['protocol'].append('pre_outcome_protocol_lock_not_documented')
    if accessed is True:
        if outcome is None or frozen is None or frozen>=outcome:
            b['protocol'].append('protocol_not_locked_before_first_outcome_access')
    elif accessed is not False:
        b['protocol'].append('outcome_access_history_unknown')
    elif p.get('first_outcome_access_at_utc') is not None:
        b['protocol'].append('contradictory_outcome_access_history')

    mapping_id=c.get('mapping_id')
    mapping=QUALIFIED_EXTERNAL_MAPPINGS.get(mapping_id) if isinstance(mapping_id,str) else None
    if not isinstance(mapping,dict):
        b['mapping'].append('no_reviewed_external_mapping_code_is_registered')
    elif (mapping.get('cohort_id')!=c.get('cohort_id') or mapping.get('release_id')!=c.get('release_id')
          or mapping.get('target_bundle_sha256')!=p.get('target_bundle_sha256')):
        b['mapping'].append('registered_mapping_does_not_match_cohort_release_and_bundle')
    blocked=any(b.values())
    result={key+'_blockers':value for key,value in b.items()}
    result.update(schema_version=SCHEMA,status='blocked_before_scoring' if blocked else 'metadata_ready_for_reviewed_mapping',
        frozen_scoring_allowed=not blocked,metadata_complete_except_mapping=not any(b[key] for key in _CATEGORIES if key!='mapping'),
        qualified_external_mapping_registry_empty=not QUALIFIED_EXTERNAL_MAPPINGS,
        legal_authorization_verified=False,source_authority_verified=False,independence_verified=False,
        chronology_verified=False,clinical_use_ready=False,person_records_processed=False,
        network_or_file_access_performed=False,
        note='Metadata assertions never certify rights, independent participants, historical assay equivalence or pre-outcome registration.')
    return result
