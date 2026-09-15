"""Check declared HRS permissions BEFORE opening any person-level file.

Rules reviewed 2026-09-15 against the HRS AI policy and August 21, 2026 DUA.
This is an operator-attestation check, not authentication or legal certification.
Never upload HRS microdata, including derived presence/absence records, to AI.
"""
from __future__ import annotations
import json

POLICY_URL = 'https://hrsdata.isr.umich.edu/data-products/ai-llm-use-policy'
DUA_URL = 'https://hrs.isr.umich.edu/sites/default/files/HRS-Sensitive-Data-Access-Use-Agreement.pdf'
PROHIBITED_CONTEXTS = {'ai_session', 'github_actions', 'public_web_app', 'third_party_service'}


def strict_json(text: str) -> dict:
    def pairs(items):
        d = {}
        for k, v in items:
            if k in d:
                raise ValueError('duplicate_metadata_key')
            d[k] = v
        return d
    def bad_constant(_):
        raise ValueError('nonfinite_metadata_value')
    result = json.loads(text, object_pairs_hook=pairs, parse_constant=bad_constant)
    if not isinstance(result, dict):
        raise ValueError('metadata_object_required')
    return result


def text(x):
    return isinstance(x, str) and bool(x.strip()) and x.strip().lower() not in {'unknown', 'pending', 'none'}


def check_access(c: dict) -> dict:
    if not isinstance(c, dict):
        raise ValueError('metadata_object_required')
    blockers = []
    kind = c.get('data_class')
    context = c.get('execution_context')
    if kind == 'public_documentation':
        if c.get('contains_person_records') is not False:
            blockers.append('documentation_must_exclude_person_records')
    elif kind == 'synthetic_fixture':
        if c.get('generated_without_hrs_microdata') is not True:
            blockers.append('synthetic_origin_not_confirmed')
    elif kind == 'reviewed_aggregate':
        for flag in ('no_person_records', 'human_disclosure_review_completed', 'approved_for_sharing'):
            if c.get(flag) is not True:
                blockers.append(flag)
        if not text(c.get('review_reference')):
            blockers.append('review_reference_required')
    elif kind == 'hrs_microdata':
        if context != 'authorized_human_local':
            blockers.append('hrs_person_records_not_for_ai_or_unattended_cloud')
        for flag in ('user_dua_approved', 'institutional_signature_confirmed',
                     'secure_environment_reviewed', 'ai_file_access_disabled',
                     'third_party_transfer_disabled', 'research_only_use',
                     'no_consulting_or_licensing_obligation', 'source_mapping_reviewed',
                     'baseline_assignment_reviewed'):
            if c.get(flag) is not True:
                blockers.append(flag)
        if not text(c.get('approval_reference')):
            blockers.append('approval_reference_required')
        requested, approved = c.get('requested_products'), c.get('approved_products')
        if not (isinstance(requested, list) and requested and all(text(s) for s in requested)
                and isinstance(approved, list) and all(text(s) for s in approved)
                and set(requested) <= set(approved)):
            blockers.append('every_requested_product_requires_approval')
        if c.get('commercial_affiliation') is True:
            if c.get('commercial_affiliation_disclosed_and_reviewed') is not True:
                blockers.append('commercial_affiliation_requires_research_plan_review')
        elif c.get('commercial_affiliation') is not False:
            blockers.append('commercial_affiliation_not_declared')
        if c.get('uses_ndi_or_cms') is True:
            # Remote access from outside the US is still outside-US access.
            if c.get('operator_country') != 'US' or c.get('execution_country') != 'US':
                blockers.append('ndi_cms_cannot_be_accessed_outside_US')
        elif c.get('uses_ndi_or_cms') is not False:
            blockers.append('outcome_linkage_product_not_declared')
        if not text(c.get('operator_country')) or not text(c.get('execution_country')):
            blockers.append('access_locations_required')
    else:
        blockers.append('unsupported_data_class')
    return {'status': 'blocked' if blockers else 'declared_conditions_met',
            'blockers': blockers, 'data_class': kind, 'policy_url': POLICY_URL,
            'legal_authorization_verified': False, 'environment_technically_verified': False,
            'clinical_use_ready': False,
            'note': 'Self-declarations do not establish valid permission or a secure environment.'}
