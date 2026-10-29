"""Preserved source/domain qualification for the separate33B cause targets.

Only aggregate provenance leaves the source adapter. Cause labels are returned
separately as in-memory arrays in the unchanged primary common-row order.
"""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd

from ..cause33 import gate


ROOT = Path(__file__).parent
PRIOR = ROOT.parent / 'cause33'
TARGET_CODES = ('002', '001')
SOURCE_FLOAT_RTOL = 1e-12
SOURCE_FLOAT_ATOL = 1e-9
SCIENCE_FIELDS = frozenset({
    'schema', 'gate', 'endpoint', 'eligible', 'common', 'per_cycle', 'causes',
    'unknown_cause', 'defined_cause_ranking', 'priority_top5', 'proceed_codes',
    'accounting', 'design', 'folds', 'weight_rule', 'fitted_models',
    'predictive_metrics_computed', 'participant_records_exported',
    'threshold_scope',
})
RUNNER_FIELDS = frozenset({
    'plan_commit', 'lock_sha256', 'runtime_python', 'prior_cohort_use',
    'comparison_models_fitted', 'independent_binary_probabilities_summed',
    'joint_integration_and_Fine_Gray', 'historical_scientific_files_verified',
    'scientific_python_sha256',
})


def _json_bytes(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':'),
                      allow_nan=False).encode('utf-8')


def _verify_wrapper(source):
    # The root-owned wrapper checks the committed33B protocol and all inherited
    # metadata/code pins before the old source parser can read participant data.
    from .retrieve import verify_wrapper
    return verify_wrapper(source)


def _load_prior(source):
    # Keep pure source/projection tests independent of the old mortality reader's
    # model-related dependencies; actual parsing uses the unchanged31 loader.
    from ..b2m31.cohort import load
    return load(source)


def _scientific_projection(value, frozen=False):
    expected = SCIENCE_FIELDS | RUNNER_FIELDS if frozen else SCIENCE_FIELDS
    if not isinstance(value, dict) or set(value) != expected:
        raise ValueError('Changed33A scientific or runner field set')
    schema = 'cause33-qualified-counts/1' if frozen else 'cause33-count-gate/1'
    if value['schema'] != schema:
        raise ValueError('Wrong33A qualification schema')
    projected = {key: value[key] for key in SCIENCE_FIELDS}
    projected['schema'] = 'cause33-count-gate/1'
    return projected


def _compare_scientific(actual, reference):
    """Keep types/counts exact; allow only fixed float reduction roundoff."""
    diagnostics = dict(floating_values_checked=0, floating_values_different=0,
                       maximum_float_abs_difference=0.)

    def compare(a, b):
        if type(a) is not type(b):
            raise ValueError('33A scientific aggregate type changed')
        if isinstance(b, dict):
            if set(a) != set(b) or not all(type(key) is str for key in b):
                raise ValueError('33A scientific aggregate fields changed')
            for key in b:
                compare(a[key], b[key])
        elif isinstance(b, list):
            if len(a) != len(b):
                raise ValueError('33A scientific aggregate length changed')
            for item, expected in zip(a, b):
                compare(item, expected)
        elif type(b) is float:
            if (not math.isfinite(a) or not math.isfinite(b)
                    or not math.isclose(a, b, rel_tol=SOURCE_FLOAT_RTOL,
                                        abs_tol=SOURCE_FLOAT_ATOL)):
                raise ValueError('33A scientific floating aggregate changed')
            diagnostics['floating_values_checked'] += 1
            diagnostics['floating_values_different'] += int(a != b)
            diagnostics['maximum_float_abs_difference'] = max(
                diagnostics['maximum_float_abs_difference'], abs(a - b))
        elif type(b) in (str, int, bool, type(None)):
            if a != b:
                raise ValueError('33A scientific count/label/fingerprint changed')
        else:
            raise ValueError('Non-JSON33A scientific aggregate')

    compare(actual, reference)
    return diagnostics


def reproduce_gate(eligible, common, source_audit, rules, definitions, frozen):
    """Requalify and compare all retained scientific aggregates before fitting."""
    result = gate.qualify(eligible, common, source_audit, rules, definitions)
    condition = result['gate']
    if (condition['source_endpoint_qualified'] is not True
            or condition['cause_labels_complete_for_common_horizon'] is not True
            or condition['later_comparisons_allowed'] is not True
            or condition['blockers'] or condition['support_blockers']):
        raise ValueError('Inherited33A source/endpoint gate blocks all comparisons')
    if result['proceed_codes'] != list(TARGET_CODES):
        raise ValueError('Changed33A supported target causes')
    current = _scientific_projection(result)
    reference = _scientific_projection(frozen, frozen=True)
    diagnostic = _compare_scientific(current, reference)
    proof = dict(exact_structure_and_counts=True,
        numeric_tolerance=dict(rtol=SOURCE_FLOAT_RTOL, atol=SOURCE_FLOAT_ATOL),
        **diagnostic,
        recomputed_scientific_sha256=hashlib.sha256(_json_bytes(current)).hexdigest(),
        frozen_scientific_sha256=hashlib.sha256(_json_bytes(reference)).hexdigest(),
        compared_fields=sorted(SCIENCE_FIELDS),
        expected_proceed_codes=list(TARGET_CODES))
    return result, proof


def _common_mask(eligible, common):
    if isinstance(common, pd.Series):
        if (not common.index.equals(eligible.index) or common.isna().any()
                or not pd.api.types.is_bool_dtype(common.dtype)):
            raise ValueError('Exact boolean common mask/index required')
        common = common.to_numpy(dtype=bool)
    else:
        common = np.asarray(common)
    if common.shape != (len(eligible),) or common.dtype != bool:
        raise ValueError('Exact boolean common mask required')
    return common


def cause_targets(eligible, common):
    """Return002/001 labels only in common-row order, never export them.

    This endpoint check remains strict if a caller changes the in-memory frame
    after source qualification. Full source/domain checks are mandatory in load.
    """
    if (not isinstance(eligible, pd.DataFrame) or not len(eligible)
            or not eligible.columns.is_unique
            or not {'dead', 'cause', 't_exam'} <= set(eligible)):
        raise ValueError('Eligible endpoint frame required')
    common = _common_mask(eligible, common)
    values = {}
    for name in ('dead', 'cause', 't_exam'):
        series = eligible[name]
        if (not pd.api.types.is_numeric_dtype(series.dtype)
                or pd.api.types.is_bool_dtype(series.dtype)):
            raise ValueError('Numeric endpoint fields required')
        values[name] = series.to_numpy(dtype=float, na_value=np.nan)
    dead, cause, time = (values[key] for key in ('dead', 'cause', 't_exam'))
    known = ~np.isnan(cause)
    if (not np.isin(dead, [0., 1.]).all()
            or not np.isfinite(time).all() or (time < 0).any()
            or not (np.isfinite(cause[known]).all()
                    and (cause[known] == np.floor(cause[known])).all()
                    and np.isin(cause[known], np.arange(1, 11)).all())
            or ((dead == 0) & known).any()
            or ((dead == 0) & (time < 96)).any()):
        raise ValueError('Invalid33B status/cause/followup')
    early = (dead == 1) & (time <= 96)
    if (common & early & ~known).any():
        raise ValueError('Unknown common early cause blocks all targets')
    return {code: (early & (cause == int(code)))[common].copy()
            for code in TARGET_CODES}


def load(source):
    """Return eligible, common, aggregate audit and source provenance.

    Participant arrays are not embedded in either exportable dictionary. This
    method never changes masks, weights, design, causes or the inherited folds.
    """
    source = Path(source)
    base, wrapper = _verify_wrapper(source)
    contract = json.loads((ROOT / 'CONTRACT.json').read_text())
    for name, expected in [('source_projection_float_rtol', SOURCE_FLOAT_RTOL),
                           ('source_projection_float_atol', SOURCE_FLOAT_ATOL)]:
        if type(contract.get(name)) is not float or contract[name] != expected:
            raise ValueError('Source projection tolerance differs from fixed33B contract')
    eligible, common, _, inherited, loaded_base = _load_prior(source)
    if loaded_base != base:
        raise ValueError('Source manifest changed between verification and parsing')
    rules = json.loads((PRIOR / 'SUPPORT_RULES.json').read_text())
    definitions = json.loads((PRIOR / 'CAUSE_DEFINITIONS.json').read_text())
    frozen_path = PRIOR / 'GATE_RESULTS.json'
    frozen_bytes = frozen_path.read_bytes()
    frozen = json.loads(frozen_bytes)
    result, proof = reproduce_gate(eligible, common, inherited, rules, definitions, frozen)
    # Exercise the strict endpoint policy before returning anything to fitting.
    # Labels remain temporary in-memory arrays and are not placed in the audit.
    cause_targets(eligible, common)
    latest_base, latest_wrapper = _verify_wrapper(source)
    if latest_base != base or latest_wrapper != wrapper or frozen_path.read_bytes() != frozen_bytes:
        raise ValueError('Protocol/source provenance changed during qualification')
    frozen_sha = hashlib.sha256(frozen_bytes).hexdigest()
    audit = dict(schema='cause33b-source-domain-audit/1',
        inherited31_source_domain_audit=inherited,
        reproduced33A_gate=result, frozen33A_gate_sha256=frozen_sha,
        scientific_gate_equality=proof, participant_records_exported=False,
        comparative_models_fitted=False, labels_exported=False)
    manifest = dict(schema='cause33b-source-manifest/1', intake=wrapper, source=base,
                    frozen33A_gate_sha256=frozen_sha)
    # Fail before any consumer can serialize a nonfinite/array-valued audit.
    _json_bytes(audit)
    _json_bytes(manifest)
    return eligible, common, audit, manifest
