"""Pure metadata qualification for a named independent replication release.

No source bytes, participant records, model fitting or access actions are read
or performed. Readiness concerns a later source-specific protocol, not approval
to download data or compare outcomes in this stage.
"""
from __future__ import annotations

from copy import deepcopy
from datetime import date
import re
from urllib.parse import urlsplit


MODES = ('algorithm_replication', 'locked_source_transport')
EVIDENCE_STATES = ('verified_true', 'verified_false', 'not_verified')
REQUIREMENTS = (
    'anonymous_person_level_release', 'processing_without_new_terms_permitted',
    'independent_programme_source', 'overlap_history_assessment',
    'population_age_domain_declared', 'same_visit_routine11',
    'cys_absolute_interpretable', 'b2m_absolute_interpretable',
    'baseline_origin_aligned', 'underlying_heart_mapping',
    'underlying_cancer_mapping', 'known_other_competing_and_vital_status',
    'full96_outcome_contract', 'sampling_design_contract', 'joint_domain_auditable',
)
TRANSPORT_REQUIREMENTS = ('frozen_source_scorer_artifact', 'transport_assay_unit_mapping')
SOURCE_KINDS = ('person_level_release', 'metadata_only', 'aggregate_publication')
ACCESS_STATES = ('public_anonymous', 'controlled_not_requested', 'provider_prohibited_ai',
                 'terms_or_registration_required', 'unknown')
SOURCE_REFERENCE_KINDS = ('official_study_metadata', 'official_access_terms',
    'official_codebook', 'primary_publication', 'programme_history', 'methodology_reference')
ACCESS_ACTION = dict(public_anonymous='not_requested', controlled_not_requested='not_requested',
    provider_prohibited_ai='prohibited', terms_or_registration_required='not_granted',
    unknown='not_verified')
ACCESS_GRANT = dict(public_anonymous='not_required_for_anonymous_public_release',
    controlled_not_requested='not_granted', provider_prohibited_ai='prohibited',
    terms_or_registration_required='not_granted', unknown='not_verified')
ACTIVATION_STATUS = 'blocked_until_named_qualified_release_and_preoutcome_mappings_precision_plan'


def _keys(value, required, label, optional=()):
    if (type(value) is not dict or not set(required) <= set(value)
            or not set(value) <= set(required) | set(optional)):
        raise ValueError('Wrong ' + label + ' field set')


def _text(value, label):
    if type(value) is not str or not value.strip():
        raise ValueError('Nonempty text required: ' + label)


def _text_list(value, label):
    if type(value) is not list:
        raise ValueError('Text list required: ' + label)
    for item in value:
        _text(item, label)


def _digest(value):
    if type(value) is not str or re.fullmatch(r'[0-9a-f]{64}', value) is None:
        raise ValueError('Lowercase SHA256 metadata digest required')


def _https(value):
    _text(value, 'source URL')
    try:
        parsed = urlsplit(value)
        parsed.port  # Validate a declared port without imposing a provider port.
        valid = (parsed.scheme == 'https' and bool(parsed.hostname)
                 and parsed.username is None and parsed.password is None
                 and not any(char.isspace() or ord(char) < 32 or ord(char) == 127
                             for char in value) and '\\' not in value)
    except ValueError:
        valid = False
    if not valid:
        raise ValueError('HTTPS source URL without credentials required')


def _sources(manifest):
    """Validate only declared metadata; never resolve a URL or open a file."""
    required = ('schema', 'checked_at_UTC', 'sources', 'prior_catalog', 'review_sha256',
                'original_source_bytes_saved', 'new_outcomes_read',
                'person_level_downloads_performed')
    _keys(manifest, required, 'source manifest')
    if manifest['schema'] != 'replication34-sources/1':
        raise ValueError('Wrong source manifest schema')
    _text(manifest['checked_at_UTC'], 'metadata review date')
    try:
        if date.fromisoformat(manifest['checked_at_UTC']).isoformat() != manifest['checked_at_UTC']:
            raise ValueError
    except ValueError:
        raise ValueError('ISO UTC metadata review date required') from None
    for field in ('original_source_bytes_saved', 'new_outcomes_read', 'person_level_downloads_performed'):
        if manifest[field] is not False:
            raise ValueError('Metadata-only source manifest cannot claim source/outcome/download actions')
    _keys(manifest['prior_catalog'], ('path', 'sha256'), 'prior catalogue')
    _text(manifest['prior_catalog']['path'], 'prior catalogue path')
    _digest(manifest['prior_catalog']['sha256'])
    reviews = manifest['review_sha256']
    if type(reviews) is not dict:
        raise ValueError('Review fingerprint dictionary required')
    for name, digest in reviews.items():
        _text(name, 'review filename'); _digest(digest)
    if type(manifest['sources']) is not list:
        raise ValueError('Source metadata list required')
    sources = {}
    for source in manifest['sources']:
        _keys(source, ('id', 'url', 'kind', 'review_file', 'observation_status',
                       'original_source_sha256'), 'source reference', optional=('detail',))
        _text(source['id'], 'source identifier')
        if source['id'] in sources:
            raise ValueError('Duplicate source identifier')
        _https(source['url'])
        if source['kind'] not in SOURCE_REFERENCE_KINDS:
            raise ValueError('Unknown source reference kind')
        _text(source['review_file'], 'source review filename')
        if source['review_file'] not in reviews:
            raise ValueError('Source review lacks declared review fingerprint')
        _text(source['observation_status'], 'source observation status')
        if source['original_source_sha256'] is not None:
            raise ValueError('Metadata ledger cannot claim an original source-byte digest')
        if 'detail' in source:
            _text(source['detail'], 'source metadata detail')
        sources[source['id']] = source
    return sources


def _candidate(candidate, sources):
    _keys(candidate, ('candidate_id', 'release_id', 'source_kind', 'access_state',
                      'checks', 'user_actions'), 'candidate')
    _text(candidate['candidate_id'], 'candidate identifier')
    _text(candidate['release_id'], 'named release/visit identifier')
    if candidate['source_kind'] not in SOURCE_KINDS:
        raise ValueError('Unknown candidate source kind')
    if candidate['access_state'] not in ACCESS_STATES:
        raise ValueError('Unknown candidate access state')
    _text_list(candidate['user_actions'], 'human next actions')
    _keys(candidate['checks'], REQUIREMENTS + TRANSPORT_REQUIREMENTS, 'candidate checks')
    for name, evidence in candidate['checks'].items():
        _keys(evidence, ('state', 'source_ids', 'reason'), 'evidence check')
        if evidence['state'] not in EVIDENCE_STATES:
            raise ValueError('Unknown evidence state')
        _text(evidence['reason'], 'evidence reason')
        _text_list(evidence['source_ids'], 'evidence source identifiers')
        if len(set(evidence['source_ids'])) != len(evidence['source_ids']):
            raise ValueError('Duplicate evidence source identifier')
        if any(identifier not in sources for identifier in evidence['source_ids']):
            raise ValueError('Unresolved evidence source identifier')
        if evidence['state'] == 'verified_true' and not evidence['source_ids']:
            raise ValueError('Verified evidence requires a resolved source reference')


def _mode(mode):
    if type(mode) is not str or mode not in MODES:
        raise ValueError('Unknown replication qualification mode')


def _qualify(candidate, sources, mode):
    _candidate(candidate, sources)
    required = REQUIREMENTS + (TRANSPORT_REQUIREMENTS if mode == 'locked_source_transport' else ())
    checks = candidate['checks']
    blockers = [dict(check=name, **deepcopy(checks[name])) for name in required
                if checks[name]['state'] != 'verified_true']
    kind, access = candidate['source_kind'], candidate['access_state']
    if kind != 'person_level_release':
        blockers.append(dict(check='source_kind', state='verified_false', source_ids=[],
            reason='Declared ' + kind + ' is not a person-level release.'))
    if access != 'public_anonymous':
        reasons = dict(controlled_not_requested='Controlled access has not been requested or granted.',
            provider_prohibited_ai='Provider prohibits the proposed AI processing.',
            terms_or_registration_required='New terms or registration are required and have not been granted.',
            unknown='Anonymous public processing access has not been verified.')
        blockers.append(dict(check='access_state',
            state='not_verified' if access == 'unknown' else 'verified_false',
            source_ids=[], reason=reasons[access]))
    ready = not blockers
    return dict(schema='replication34-candidate-qualification/1',
        candidate_id=candidate['candidate_id'], release_id=candidate['release_id'], mode=mode,
        source_kind=kind, access_state=access, source_ready=ready,
        source_ready_for_protocol=ready, required_checks=list(required),
        checks=deepcopy(checks), blockers=blockers, user_actions=list(candidate['user_actions']),
        access_action_state=ACCESS_ACTION[access], access_grant_state=ACCESS_GRANT[access],
        access_actions_performed=False, person_level_download_authorized=False,
        outcome_comparison_authorized=False, protocol_activation_status=ACTIVATION_STATUS,
        precision_adequacy='not_established',
        algorithm_precision_unqualified=(mode == 'algorithm_replication'),
        qualification_scope='metadata_readiness_only_no_download_fit_or_activation')


def qualify(candidate, source_manifest, mode='algorithm_replication'):
    """Qualify declared, cited metadata; never grant access or activate a study."""
    _mode(mode)
    return _qualify(candidate, _sources(source_manifest), mode)


def qualify_all(candidates, source_manifest, mode='algorithm_replication'):
    """Return complete aggregate qualifications in the supplied candidate order."""
    _mode(mode)
    sources = _sources(source_manifest)
    if type(candidates) is not list:
        raise ValueError('Candidate list required')
    identifiers = set()
    for candidate in candidates:
        _candidate(candidate, sources)
        identifier = candidate['candidate_id']
        if identifier in identifiers:
            raise ValueError('Duplicate candidate identifier')
        identifiers.add(identifier)
    results = [_qualify(candidate, sources, mode) for candidate in candidates]
    ready = [result['candidate_id'] for result in results if result['source_ready']]
    return dict(schema='replication34-source-qualifications/1', mode=mode,
        qualifications=results, candidate_count=len(results),
        source_ready_candidate_ids=ready, source_ready=bool(ready),
        qualified_candidate_count=len(ready), precision_adequacy='not_established',
        algorithm_precision_unqualified=(mode == 'algorithm_replication'),
        outcome_comparison_authorized=False, person_level_download_authorized=False,
        access_actions_performed=False)
