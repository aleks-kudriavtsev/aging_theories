"""Strict in-memory research scorer; no participant reader or scoring CLI.

Only the documented training-release representation is initially supported.
Passing a schema does not authenticate provenance, authorize access, or validate
an assay. Arrays returned to an authorized caller must not be published as
participant-level evidence. No inference-time fitting or unit conversion occurs.
"""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd

from ..binary28.cohort import engineer as clinical_engineer
from ..routine29.features import engineer, CORE
from ..routine29.model import CORE_NAMES, PROTEIN_NAMES, transform_features, predict_ridge

ROOT = Path(__file__).parent
ARMS = ('routine', 'routine_all3')
ARTIFACT_FIELDS = {'schema', 'target', 'mapping_id', 'contract_sha256', 'plan_sha256',
    'lock_sha256', 'plan_commit', 'model_family', 'fit_scope', 'selection_scope',
    'models', 'training', 'source_manifest_sha256', 'source_content_sha256',
    'source_fingerprints', 'source_policy', 'scientific_python_sha256', 'versions',
    'clinical_use_ready', 'independent_external_validation', 'participant_records_exported'}


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _pairs(items):
    result = {}
    for key, value in items:
        if key in result:
            raise ValueError('Duplicate JSON key')
        result[key] = value
    return result


def strict_json(path):
    path = Path(path)
    if path.stat().st_size > 2_000_000:
        raise ValueError('Metadata size limit exceeded')
    def finite_float(token):
        value = float(token)
        if not math.isfinite(value):
            raise ValueError('Nonfinite JSON number')
        return value
    return json.loads(path.read_text(), object_pairs_hook=_pairs, parse_float=finite_float,
                      parse_constant=lambda _: (_ for _ in ()).throw(ValueError('Nonfinite JSON')))


def check_study():
    lock = strict_json(ROOT / 'LOCK.json')
    for name, expected in lock['metadata_sha256'].items():
        path = ROOT.parent / name if '/' in name else ROOT / name
        if digest(path) != expected:
            raise ValueError('Frozen study metadata changed')
    return strict_json(ROOT / 'CONTRACT.json')


def _number(value):
    if type(value) not in (float, int):
        return False
    try:
        return math.isfinite(value)
    except (ValueError, OverflowError):
        return False


def _vector(value, n, positive=False):
    if (not isinstance(value, list) or len(value) != n
            or any(not _number(v) for v in value)):
        raise ValueError('Invalid frozen numeric vector')
    a = np.asarray(value, dtype=float)
    if positive and (a <= 0).any():
        raise ValueError('Frozen scales must be positive')
    return a


def validate_artifact(artifact):
    contract = check_study()
    if (not isinstance(artifact, dict) or set(artifact) - ARTIFACT_FIELDS
            or artifact.get('schema') != contract['artifact_schema']
            or artifact.get('target') != contract['target']
            or artifact.get('mapping_id') != contract['mapping_id']
            or artifact.get('contract_sha256') != digest(ROOT / 'CONTRACT.json')
            or artifact.get('plan_sha256') != digest(ROOT / 'FREEZE_PLAN.md')
            or not isinstance(artifact.get('models'), dict)
            or set(artifact['models']) != set(ARMS)):
        raise ValueError('Frozen artifact contract mismatch')
    for arm in ARMS:
        row = artifact['models'][arm]
        if not isinstance(row, dict) or not isinstance(row.get('preprocessing'), dict) or not isinstance(row.get('fit'), dict):
            raise ValueError('Invalid frozen model structure')
        if set(row) - {'preprocessing', 'fit', 'feature_names', 'inner_brier_by_penalty',
                       'selected_penalty', 'inner_folds_aggregate', 'inner_fit_diagnostics'}:
            raise ValueError('Unknown frozen model fields')
        spec, fit = row['preprocessing'], row['fit']
        if (set(spec) != {'model_name', 'core_indices', 'protein_transforms', 'mean', 'scale', 'constant'}
                or set(fit) - {'intercept', 'coefficients', 'gradient_max', 'objective',
                               'iterations', 'penalty', 'convergence_status', 'normalized_weight_sum'}):
            raise ValueError('Unknown preprocessing or fit fields')
        n = 11 if arm == 'routine' else 14
        transforms = [] if arm == 'routine' else [dict(column=i, kind='identity_log2') for i in range(3)]
        names = list(CORE_NAMES) + ([] if arm == 'routine' else list(PROTEIN_NAMES))
        if (spec.get('model_name') != arm or spec.get('core_indices') != list(range(11))
                or any(type(v) is not int for v in spec['core_indices'])
                or spec.get('protein_transforms') != transforms
                or any(type(t.get('column')) is not int for t in spec['protein_transforms'])
                or row.get('feature_names') != names):
            raise ValueError('Unapproved frozen feature transformation')
        _vector(spec.get('mean'), n)
        _vector(spec.get('scale'), n, positive=True)
        if (not isinstance(spec.get('constant'), list) or len(spec['constant']) != n
                or any(type(v) is not bool for v in spec['constant'])):
            raise ValueError('Invalid constant-feature flags')
        _vector(fit.get('coefficients'), n)
        if (not _number(fit.get('intercept')) or type(fit.get('penalty')) not in (float, int)
                or fit['penalty'] not in (1., 10., 100.)
                or not _number(row.get('selected_penalty'))
                or row.get('selected_penalty') != fit['penalty']):
            raise ValueError('Invalid frozen fit')
    return artifact


def load_artifact(path, manifest_path):
    """Check exact artifact bytes against a separately pinned freeze manifest.

    Hashes detect changes relative to that manifest, not malicious replacement
    of both files. Repository commits establish the reviewed reference.
    """
    manifest = strict_json(manifest_path)
    if (manifest.get('schema') != 'freeze30-manifest/1'
            or manifest.get('artifact_sha256') != digest(path)
            or manifest.get('contract_sha256') != digest(ROOT / 'CONTRACT.json')
            or manifest.get('plan_sha256') != digest(ROOT / 'FREEZE_PLAN.md')):
        raise ValueError('Freeze manifest mismatch')
    return validate_artifact(strict_json(path))


def input_features(document):
    """Validate canonical released measurements before engineering features."""
    contract = check_study()
    keys = {'schema', 'contract_sha256', 'mapping_id', 'units', 'records'}
    if not isinstance(document, dict) or set(document) != keys:
        raise ValueError('Exact scoring document schema required')
    if (document['schema'] != contract['input_schema']
            or document['contract_sha256'] != digest(ROOT / 'CONTRACT.json')
            or document['mapping_id'] != contract['mapping_id']
            or document['units'] != contract['units']):
        raise ValueError('Unapproved schema, units or measurement mapping')
    rows = document['records']
    if not isinstance(rows, list) or not 1 <= len(rows) <= 10_000:
        raise ValueError('Nonempty bounded scoring batch required')
    fields = contract['fields']
    for row in rows:
        if not isinstance(row, dict) or set(row) != set(fields):
            raise ValueError('Exact measurement fields required')
        for name, spec in fields.items():
            value = row[name]
            if name in ('sex', 'smoking'):
                if not isinstance(value, str) or value not in spec['values']:
                    raise ValueError('Unknown measurement category')
                continue
            if (not _number(value) or (spec.get('type') == 'integer' and type(value) is not int)
                    or ('values' in spec and value not in spec['values'])
                    or ('minimum' in spec and value < spec['minimum'])
                    or ('maximum' in spec and value > spec['maximum'])
                    or value in spec.get('sentinels', [])):
                raise ValueError('Invalid or unsupported measurement')
            if 'values' not in spec and value <= 0:
                raise ValueError('Positive reported measurement required')
    d = pd.DataFrame(rows)
    raw = pd.DataFrame(dict(
        HSAGEIR=d.age_years, HSSEX=d.sex.map({'male': 1., 'female': 2.}),
        HAR1=np.where(d.smoking.eq('never'), 2., 1.),
        HAR3=np.where(d.smoking.eq('current'), 1., 2.),
        HAZMNK1R=d.sbp_mmhg, HAZNOK1R=d.sbp_measurement_count,
        CEP=d.creatinine_historical_mg_dl, GHP=d.hba1c_percent, GHPMETH=d.hba1c_method,
        URP=d.urine_creatinine_mg_dl, UBP=d.urine_albumin_ug_ml,
        BMPBMI=d.bmi_kg_m2, BMPWTFLG=d.bmi_weight_flag, BMPHTFLG=d.bmi_height_flag))
    try:
        engineered = engineer(clinical_engineer(raw))
        core = engineered[list(CORE)].to_numpy(float)
        proteins = np.log2(d[['cystatin_c_reported_mg_l', 'beta2_microglobulin_reported_mg_l',
                              'beta_trace_protein_reported_mg_l']].to_numpy(float))
        if not np.isfinite(core).all() or not np.isfinite(proteins).all():
            raise ValueError('Invalid feature domain')
    except (ValueError, ArithmeticError, OverflowError):
        raise ValueError('Measurements cannot form the frozen feature domain') from None
    return core, proteins


def score(document, artifact):
    """Return paired probabilities in memory; never write or log participant data."""
    validate_artifact(artifact)
    core, proteins = input_features(document)
    return {arm: predict_ridge(transform_features(core, proteins, artifact['models'][arm]['preprocessing']),
                              artifact['models'][arm]['fit']) for arm in ARMS}
