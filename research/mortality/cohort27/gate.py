"""Metadata-only source/endpoint/assay gate for NHANES III qualification.

No file is opened and no legal right is established by these checks. A passing
qualification gate permits only numeric-availability/endpoint counting, never
model validation, quantification-limit adjudication or clinical deployment.
"""
from __future__ import annotations

import re
from urllib.parse import urlsplit

COHORT = 'NHANESIII_FIRST_EXAM'
PUBLIC_HOSTS = frozenset({'www.cdc.gov', 'wwwn.cdc.gov', 'ftp.cdc.gov'})
CONTEXTS = frozenset({'ai_session', 'github_actions', 'authorized_human_local'})
ASSAYS = {
    'CST3': {'variable': 'SSCYSTAT', 'uniprot': 'P01034'},
    'B2M': {'variable': 'SSB2M', 'uniprot': 'P61769'},
    'BTP': {'variable': 'SSBTP', 'uniprot': 'P41222'},
}
AVAILABILITY_POLICY = 'released_values_numeric_availability_only_not_LOQ_adjudicated'
MODEL_FLAGS = (
    'required_features_mapping_verified',
    'clinical_semantics_verified',
    'assay_transport_reviewed',
)


def _text(value):
    return (isinstance(value, str) and bool(value.strip())
            and value.strip().lower() not in {'unknown', 'pending', 'none', 'not_reported'})


def _official_url(value):
    if not _text(value):
        return False
    try:
        parsed = urlsplit(value)
        return (parsed.scheme == 'https' and parsed.hostname in PUBLIC_HOSTS
                and parsed.username is None and parsed.password is None
                and parsed.port in (None, 443) and bool(parsed.path)
                and not parsed.fragment)
    except (ValueError, TypeError):
        return False


def _mapping(contract, key):
    value = contract.get(key, {})
    if not isinstance(value, dict):
        raise ValueError(key + '_must_be_metadata_object')
    return value


def audit(contract):
    """Assess declared metadata before any participant-file read.

    Self-reported source review and official documentation locators do not
    authenticate authorization or establish measurement interchangeability.
    Only the known CDC first-examination release is supported for real intake.
    The separate synthetic flag must never be interpreted as empirical access.
    """
    if not isinstance(contract, dict):
        raise ValueError('metadata_object_required')
    blockers = []
    unresolved = []
    cohort = contract.get('cohort_id')
    kind = contract.get('data_class')
    context = contract.get('execution_context')
    synthetic = kind == 'synthetic_fixture'

    if not isinstance(context, str) or context not in CONTEXTS:
        blockers.append('unsupported_execution_context')
    if (isinstance(cohort, str) and cohort.upper().startswith('HRS')
            and isinstance(context, str) and context in {'ai_session', 'github_actions'}):
        blockers.append('hrs_person_records_prohibited_in_ai_or_unattended_cloud')
    if cohort != COHORT:
        blockers.append('unsupported_cohort_for_this_importer')
    if kind == 'controlled':
        # This public-only adapter does not upgrade a permission attestation
        # into permission to process controlled microdata in an AI session.
        blockers.append('controlled_data_not_supported_by_public_adapter')
        if contract.get('access_approved') is not True:
            blockers.append('controlled_access_not_approved')
    elif kind == 'synthetic_fixture':
        if contract.get('generated_without_microdata') is not True:
            blockers.append('synthetic_origin_not_confirmed')
    elif kind != 'public_unrestricted':
        blockers.append('unsupported_data_class')
    if contract.get('source_reviewed') is not True:
        blockers.append('public_source_review_not_confirmed')
    if not _official_url(contract.get('source_permission_evidence')):
        blockers.append('official_public_source_evidence_required')
    if contract.get('baseline_visit') != 'first_exam':
        blockers.append('first_exam_required_second_exam_not_interchangeable')

    endpoint = _mapping(contract, 'endpoint')
    for key, value in [('kind', 'all_cause_death'), ('baseline_origin', 'exam'),
                       ('followup_unit', 'months')]:
        if endpoint.get(key) != value:
            blockers.append('endpoint_' + key + '_unresolved')
    for key in ('followup_after_baseline_verified', 'linkage_same_person_verified'):
        if endpoint.get(key) is not True:
            blockers.append('endpoint_' + key)
    horizon = endpoint.get('horizon_years')
    if type(horizon) not in (int, float) or horizon != 8:
        blockers.append('locked_eight_year_qualification_horizon_required')

    sampling = _mapping(contract, 'sampling')
    if sampling.get('weight_variable') != 'WTCYPEX6':
        blockers.append('first_exam_protein_subsample_weight_WTCYPEX6_required')
    if sampling.get('subsample_weights_documented') is not True:
        blockers.append('protein_subsample_weights_not_documented')

    assays = _mapping(contract, 'assays')
    for key, identity in ASSAYS.items():
        row = assays.get(key)
        if not isinstance(row, dict):
            blockers.append(key + ':assay_metadata_missing')
            continue
        for field, value in dict(identity, unit='mg/L', matrix='serum',
                                 scale='absolute_concentration', baseline_visit='first_exam').items():
            if row.get(field) != value:
                blockers.append(key + ':' + field + '_mismatch')
        if not _official_url(row.get('source_locator')):
            blockers.append(key + ':official_assay_documentation_required')
        if not _text(row.get('method')):
            blockers.append(key + ':assay_method_unresolved')
        if row.get('quantification_policy') != AVAILABILITY_POLICY:
            blockers.append(key + ':numeric_availability_policy_required')

    model = _mapping(contract, 'frozen_model')
    for flag in MODEL_FLAGS:
        if model.get(flag) is not True:
            unresolved.append('frozen_model:' + flag)
    digest = model.get('baseline_model_sha256')
    if not isinstance(digest, str) or re.fullmatch('[0-9a-f]{64}', digest) is None:
        unresolved.append('frozen_model:valid_model_sha256_required')
    # Even complete operator assertions cannot authorize a new evaluation in a
    # qualification-only protocol; a separately frozen evaluation plan is needed.
    unresolved.append('frozen_model:separate_evaluation_protocol_required')
    unresolved.append('assay:quantification_and_transport_not_established_by_numeric_availability')

    ready = not blockers
    return {
        'status': 'qualification_metadata_ready' if ready else 'blocked_before_participant_read',
        'qualification_import_allowed': ready and not synthetic,
        'synthetic_fixture_import_allowed': ready and synthetic,
        'frozen_model_evaluation_allowed': False,
        'blockers': blockers,
        'unresolved_metadata': unresolved,
        'qualification_scope': 'numeric_availability_and_endpoint_counts_only',
        'quantified_or_QC_usable_measurements_established': False,
        'notes': [
            'First-exam protein codebook field units are mg/L; BTP prose contains a conflicting assay-range mg/dL statement.',
            'No automatic range clipping, mg/dL conversion, LOQ classification or assay harmonization is authorized.',
            'Repeat-examination convenience samples and their absent survey weights cannot replace first examination.',
            'A successful metadata gate is not independent validation or evidence of clinical readiness.',
        ],
        'legal_authorization_verified': False,
        'clinical_use_ready': False,
        'participant_records_processed': False,
    }


def require_qualification(contract):
    """Raise before opening real files unless the public intake gate passes."""
    result = audit(contract)
    if result['qualification_import_allowed'] is not True:
        raise ValueError('qualification_import_blocked: ' + '; '.join(result['blockers']))
    return result
