"""Pure source/count qualification, without fitting or predictive metrics.

Inputs are an already source-verified b2m31 eligible frame and its outcome-free
primary mask. The returned object contains aggregate counts and public-design
fingerprints only, never participant keys, individual labels or predictions.
"""
from __future__ import annotations

import hashlib
import json

import numpy as np
import pandas as pd


CORE = ('age_linear', 'age_quadratic', 'male', 'current', 'former', 'sbp',
        'bmi', 'hba1c', 'egfr_low', 'egfr_high', 'log_uacr')
CODES = tuple(f'{number:03d}' for number in range(1, 11))
THRESHOLDS = ('event_rows_min', 'noncase_rows_min', 'event_psus_min',
              'outer_evaluation_event_rows_min', 'outer_evaluation_noncase_rows_min',
              'outer_training_event_rows_min', 'outer_training_noncase_rows_min',
              'inner_evaluation_event_rows_min', 'inner_evaluation_noncase_rows_min',
              'inner_training_event_rows_min', 'inner_training_noncase_rows_min')
REQUIRED = ('SEQN', 'cycle', 'SDDSRVYR', 'RIDAGEYR', 'RIAGENDR', 'RIDSTATR',
            'elig', 'dead', 'cause', 't_exam', 'SDMVSTRA', 'SDMVPSU',
            'WTSSCB4Y', 'WTSSCB2Y', 'weight', 'SSCYST', 'SSB2M',
            'SSCYSTL', 'SSB2ML') + CORE


def _numeric(frame, name):
    series = frame[name]
    if (not pd.api.types.is_numeric_dtype(series.dtype)
            or pd.api.types.is_bool_dtype(series.dtype)):
        raise ValueError('Nonnumeric or boolean source field: ' + name)
    try:
        return series.to_numpy(dtype=float, na_value=np.nan)
    except (TypeError, ValueError):
        raise ValueError('Invalid numeric source field: ' + name) from None


def _integer(values, name):
    if not (np.isfinite(values).all() and np.equal(values, np.floor(values)).all()):
        raise ValueError('Missing, nonfinite or fractional source field: ' + name)


def _inventory(strata, psu, mask):
    return {str(int(s)): sorted({int(p) for p in psu[mask & (strata == s)]})
            for s in sorted(set(strata[mask]))}


def _partition_map(keys, salt, folds):
    """Hash-sort public masked clusters, exactly matching the retained31 maps."""
    ordered = sorted(set(keys), key=lambda key: hashlib.sha256((salt + ':' + key).encode()).hexdigest())
    if len(ordered) < folds:
        raise ValueError('Insufficient masked clusters for fixed partition')
    return {key: index % folds for index, key in enumerate(ordered)}


def _fingerprint(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def _contracts(rules, definitions):
    if (not isinstance(rules, dict) or not isinstance(definitions, dict)
            or definitions.get('schema') != 'cause33-definitions/1'
            or definitions.get('underlying_cause_field') != 'UCOD_LEADING'
            or definitions.get('not_contributing_flags') is not True
            or definitions.get('supported_baseline_start_years') != [1999, 2001, 2003]
            or definitions.get('horizon_months') != 96
            or definitions.get('unknown_early_cause_policy') != 'block_all_comparative_fits_not_all_negative'
            or set(definitions.get('groups', {})) != set(CODES)):
        raise ValueError('Invalid fixed public cause definitions')
    for code in CODES:
        group = definitions['groups'][code]
        if (not isinstance(group, dict)
                or group.get('eligible_for_defined_clinical_ranking') is not (code != '010')):
            raise ValueError('Changed defined-cause ranking contract')
    if (rules.get('schema') != 'cause33-support/1' or rules.get('rank_maximum') != 5
            or rules.get('rank_before_support_not_replace_failed_top5') is not True
            or rules.get('rank_order') != ['weighted_events8_descending', 'raw_events8_descending', 'code_ascending']
            or rules.get('outer_folds') != 5 or rules.get('inner_folds') != 3
            or rules.get('outer_salt') != 'b2m31:outer'
            or rules.get('inner_salt_pattern') != 'b2m31:inner:{outer}'
            or rules.get('comparison_fitting_in_this_stage') is not False):
        raise ValueError('Changed fixed support or fold contract')
    for key in THRESHOLDS:
        if type(rules.get(key)) is not int or rules[key] < 0:
            raise ValueError('Invalid support threshold: ' + key)
    for key in ('known_prior_common_n', 'known_prior_common_events8'):
        if rules.get(key) is not None and (type(rules[key]) is not int or rules[key] < 0):
            raise ValueError('Invalid prior-domain count contract')


def qualify(eligible, common, source_audit, rules, definitions):
    """Validate the unchanged frame and return only source/count aggregates.

    Unknown common early causes produce a blocked report, never all-negative
    labels. Structural contract errors raise ValueError without row diagnostics.
    Optional known_prior_common_* rules are omitted only by synthetic fixtures.
    Actual source verification remains mandatory in the calling wrapper.
    """
    _contracts(rules, definitions)
    if (not isinstance(eligible, pd.DataFrame) or not len(eligible)
            or not eligible.columns.is_unique or not set(REQUIRED) <= set(eligible.columns)
            or not {'stratum', 'cluster'} <= set(eligible.columns)
            or not isinstance(source_audit, dict)):
        raise ValueError('Eligible source frame/schema required')
    if isinstance(common, pd.Series):
        if (not common.index.equals(eligible.index) or common.isna().any()
                or not pd.api.types.is_bool_dtype(common.dtype)):
            raise ValueError('Boolean common mask/index required')
        common = common.to_numpy(dtype=bool)
    else:
        common = np.asarray(common)
    if common.shape != (len(eligible),) or common.dtype != bool:
        raise ValueError('Exact one-dimensional boolean common mask required')
    common = common.copy()
    x = {name: _numeric(eligible, name) for name in set(REQUIRED)}
    for name in ('SEQN', 'cycle', 'SDDSRVYR', 'RIDAGEYR', 'RIAGENDR',
                 'RIDSTATR', 'elig', 'dead', 'SDMVSTRA', 'SDMVPSU'):
        _integer(x[name], name)
    if not ((x['SEQN'] > 0).all() and len(np.unique(x['SEQN'])) == len(eligible)):
        raise ValueError('Participant keys must be positive and unique')
    if (not np.isin(x['cycle'], [1999, 2001, 2003]).all()
            or not np.array_equal(x['SDDSRVYR'], np.where(x['cycle'] == 1999, 1, np.where(x['cycle'] == 2001, 2, 3)))):
        raise ValueError('Unsupported/coarsened baseline cycle or release mapping')
    if (not ((x['RIDAGEYR'] >= 40) & (x['RIDAGEYR'] <= 75)).all()
            or not np.isin(x['RIAGENDR'], [1, 2]).all()
            or not (x['RIDSTATR'] == 2).all() or not (x['elig'] == 1).all()
            or not np.isin(x['dead'], [0, 1]).all()):
        raise ValueError('Changed source eligibility/status domain')
    if not (np.isfinite(x['t_exam']).all() and (x['t_exam'] >= 0).all()):
        raise ValueError('Nonfinite or negative exam followup')
    if ((x['dead'] == 0) & (x['t_exam'] < 96)).any():
        raise ValueError('Survivor censored before fixed horizon')
    cause = x['cause']; missing = np.isnan(cause); known = ~missing
    if (not (np.isfinite(cause[known]).all()
             and np.equal(cause[known], np.floor(cause[known])).all()
             and np.isin(cause[known], np.arange(1, 11)).all())
            or ((x['dead'] == 0) & known).any()):
        raise ValueError('Undeclared/fractional cause or cause recorded for alive participant')
    for cycle, low, high, allowed_psu in [(1999, 1, 13, [1, 2, 3]), (2001, 14, 28, [1, 2]), (2003, 29, 43, [1, 2])]:
        mask = x['cycle'] == cycle
        if not (((x['SDMVSTRA'][mask] >= low) & (x['SDMVSTRA'][mask] <= high)).all()
                and np.isin(x['SDMVPSU'][mask], allowed_psu).all()):
            raise ValueError('Cycle and masked design mismatch')
    clusters = np.array([f'{int(s)}:{int(p)}' for s, p in zip(x['SDMVSTRA'], x['SDMVPSU'])])
    strata_names = np.array([f'NHANES:{int(s)}' for s in x['SDMVSTRA']])
    if (not np.array_equal(eligible.cluster.to_numpy(), clusters)
            or not np.array_equal(eligible.stratum.to_numpy(), strata_names)):
        raise ValueError('Noncanonical masked cluster/stratum identities')
    early_cycle = np.isin(x['cycle'], [1999, 2001])
    if (not (np.isfinite(x['WTSSCB2Y']).all() and (x['WTSSCB2Y'] > 0).all())
            or not (np.isfinite(x['WTSSCB4Y'][early_cycle]).all() and (x['WTSSCB4Y'][early_cycle] > 0).all())
            or not np.isnan(x['WTSSCB4Y'][~early_cycle]).all()):
        raise ValueError('Invalid retained native subsample weights')
    w = x['weight']
    expected_w = np.where(early_cycle, (2 / 3) * x['WTSSCB4Y'], (1 / 3) * x['WTSSCB2Y'])
    if (not (np.isfinite(w).all() and (w > 0).all()) or not np.array_equal(w, expected_w)):
        raise ValueError('Pooled weight differs from exact declared six-year formula')
    if 'weight6y' in eligible and not np.array_equal(_numeric(eligible, 'weight6y'), expected_w):
        raise ValueError('Retained six-year weight changed')
    for measurement, flag in [('SSCYST', 'SSCYSTL'), ('SSB2M', 'SSB2ML')]:
        measured = ~np.isnan(x[measurement]); observed_flag = ~np.isnan(x[flag])
        if (not np.array_equal(measured, observed_flag)
                or not (np.isfinite(x[measurement][measured]).all() and (x[measurement][measured] > 0).all())
                or not np.isin(x[flag][observed_flag], [0, 1, 2]).all()):
            raise ValueError('Invalid reported protein/detection field')
    expected_common = (np.isfinite(np.column_stack([x[name] for name in CORE])).all(axis=1)
                       & np.isfinite(x['SSCYST']) & (x['SSCYST'] > 0)
                       & np.isfinite(x['SSB2M']) & (x['SSB2M'] > 0)
                       & (x['SSCYSTL'] == 0) & (x['SSB2ML'] == 0))
    if not np.array_equal(common, expected_common):
        raise ValueError('Common mask changed from outcome-free primary flag00 domain')
    full = {str(s): [1, 2, 3] if s == 1 else [1, 2] for s in range(1, 44)}
    if (source_audit.get('full_DEMO_masked_design') != full
            or _inventory(x['SDMVSTRA'], x['SDMVPSU'], np.ones(len(w), dtype=bool)) != full
            or _inventory(x['SDMVSTRA'], x['SDMVPSU'], common) != full):
        raise ValueError('Full/eligible/common masked inventory not preserved')
    if 'eligible_masked_design' in source_audit and source_audit['eligible_masked_design'] != full:
        raise ValueError('Inherited eligible masked inventory changed')
    died = x['dead'] == 1; death8 = died & (x['t_exam'] <= 96)
    unknown8 = death8 & missing; unknown_later = died & (x['t_exam'] > 96) & missing
    for key, actual in [('known_prior_common_n', int(common.sum())),
                        ('known_prior_common_events8', int((common & death8).sum()))]:
        if rules.get(key) is not None and rules[key] != actual:
            raise ValueError('Prior common domain/event count changed')
    if not (np.isfinite(w.sum()) and np.isfinite(w @ w)):
        raise ValueError('Nonfinite survey-weight aggregates')

    outer_map = _partition_map(clusters, rules['outer_salt'], 5)
    outer = np.array([outer_map[key] for key in clusters], dtype=int)
    cells = []
    fold_diagnostics = []
    for fold in range(5):
        training = outer != fold; evaluation = outer == fold
        inner_map = _partition_map(clusters[training], rules['inner_salt_pattern'].format(outer=fold), 3)
        inner = np.array([inner_map.get(key, -1) for key in clusters], dtype=int)
        cells.extend([(dict(outer_fold=fold, level='outer', subset='training'), common & training),
                      (dict(outer_fold=fold, level='outer', subset='evaluation'), common & evaluation)])
        for validation in range(3):
            cells.extend([(dict(outer_fold=fold, inner_fold=validation, level='inner', subset='training'), common & training & (inner != validation)),
                          (dict(outer_fold=fold, inner_fold=validation, level='inner', subset='evaluation'), common & training & (inner == validation))])
        fold_diagnostics.append(dict(outer_fold=fold, training_PSU_n=len(inner_map),
            evaluation_PSU_n=len(set(clusters[evaluation])), PSU_overlap_n=0,
            inner_mapping_sha256=_fingerprint(inner_map),
            inner_evaluation_PSU_counts={str(i): sum(value == i for value in inner_map.values()) for i in range(3)}))

    def summary(mask, event, clinical_target=True):
        event = mask & event
        weights = w[event]; event_sum = float(weights.sum()); squared = float(weights @ weights)
        denominator = float(w[mask].sum()); all_death_weight = float(w[mask & death8].sum())
        event_n = int(event.sum()); unknown_n = int((mask & unknown8).sum())
        return dict(n=int(mask.sum()), events8=event_n, positive_weight_event_rows=event_n,
            noncase_rows=(int(mask.sum()) - event_n - unknown_n) if clinical_target else None,
            label_unknown_rows=unknown_n if clinical_target else None, event_PSU_n=len(set(clusters[event])),
            event_strata_n=len(set(x['SDMVSTRA'][event])),
            weighted_events8=event_sum, event_weight_squared_sum=squared,
            event_Kish_n=(event_sum * event_sum / squared) if event_n else None,
            max_event_weight8=float(weights.max()) if event_n else None,
            maximum_event_weight_fraction=float(weights.max() / event_sum) if event_n else None,
            common_weight_sum=denominator,
            weighted_incidence8=event_sum / denominator if denominator else None,
            weighted_share_of_all_deaths8=event_sum / all_death_weight if all_death_weight else None,
            raw_share_of_all_deaths8=event_n / int((mask & death8).sum()) if (mask & death8).any() else None)

    def population(mask):
        per_code = {code: int((mask & died & (cause == int(code))).sum()) for code in CODES}
        return dict(n=int(mask.sum()), all_cause_events8=int((mask & death8).sum()),
            weighted_all_cause_events8=float(w[mask & death8].sum()), weight_sum=float(w[mask].sum()),
            alive_n=int((mask & ~died).sum()), all_time_deaths_n=int((mask & died).sum()),
            all_time_cause_counts={**per_code, 'unknown': int((mask & died & missing).sum())},
            unknown_cause_deaths8=int((mask & unknown8).sum()),
            unknown_cause_deaths_after96=int((mask & unknown_later).sum()))

    unknown_common = int((common & unknown8).sum())
    blockers = ['unknown_common_cause_at_or_before96'] if unknown_common else []
    groups = {}
    for code in CODES:
        event = death8 & (cause == int(code))
        overall = summary(common, event)
        fold_cells = [dict(**where, **summary(mask, event)) for where, mask in cells]
        failures = []
        checks = [('event_rows_min', overall['positive_weight_event_rows'], {}),
                  ('noncase_rows_min', overall['noncase_rows'], {}),
                  ('event_psus_min', overall['event_PSU_n'], {})]
        for cell in fold_cells:
            where = {key: cell[key] for key in ('outer_fold', 'inner_fold', 'level', 'subset') if key in cell}
            prefix = cell['level'] + '_' + cell['subset']
            checks.extend([(prefix + '_event_rows_min', cell['positive_weight_event_rows'], where),
                           (prefix + '_noncase_rows_min', cell['noncase_rows'], where)])
        for key, observed, where in checks:
            if observed < rules[key]:
                failures.append(dict(rule=key, required=rules[key], observed=observed, **where))
        groups[code] = dict(code=code, cause_id=definitions['groups'][code]['id'], overall=overall,
            per_cycle={str(c): summary(common & (x['cycle'] == c), event) for c in (1999, 2001, 2003)},
            fold_cells=fold_cells, support=dict(threshold_check_count=len(checks),
                thresholds_passed=not failures, failures=failures),
            eligible_for_defined_clinical_ranking=(code != '010'))
    ranking = sorted(CODES[:-1], key=lambda code: (-groups[code]['overall']['weighted_events8'],
        -groups[code]['overall']['events8'], code))
    priority = ranking[:5]
    for index, code in enumerate(ranking, 1):
        groups[code]['defined_cause_rank'] = index
    groups['010']['defined_cause_rank'] = None
    for code in CODES:
        groups[code]['priority_top5'] = code in priority
        groups[code]['proceed_to_later_protocol'] = (code in priority
            and groups[code]['support']['thresholds_passed'] and not blockers)
    unknown = dict(**summary(common, unknown8, clinical_target=False),
        per_cycle={str(c): summary(common & (x['cycle'] == c), unknown8, clinical_target=False)
                   for c in (1999, 2001, 2003)},
        clinically_ranked=False, support_evaluated=False)
    accounted_n = sum(group['overall']['events8'] for group in groups.values()) + unknown['events8']
    accounted_weight = sum(group['overall']['weighted_events8'] for group in groups.values()) + unknown['weighted_events8']
    if (accounted_n != int((common & death8).sum()) or not np.isclose(accounted_weight,
            w[common & death8].sum(), rtol=1e-12, atol=1e-9)):
        raise ArithmeticError('Disjoint cause accounting failed')
    proceed_codes = [code for code in priority if groups[code]['proceed_to_later_protocol']]
    support_blockers = ['no_supported_priority_cause'] if not blockers and not proceed_codes else []
    return dict(schema='cause33-count-gate/1', gate=dict(source_structure_valid=True,
            cause_labels_complete_for_common_horizon=not blockers,
            source_endpoint_qualified=not blockers,
            later_comparisons_allowed=bool(proceed_codes), blockers=blockers,
            support_blockers=support_blockers,
            no_comparative_fits_performed=True),
        endpoint=dict(horizon_months=96, underlying_cause_field='UCOD_LEADING',
            competing_deaths_are_known_noncase=True, unknown_early_deaths_are_not_noncase=True),
        eligible=population(np.ones(len(w), dtype=bool)), common=population(common),
        per_cycle={str(c): dict(eligible=population(x['cycle'] == c),
            common=population(common & (x['cycle'] == c))) for c in (1999, 2001, 2003)},
        causes=groups, unknown_cause=unknown,
        defined_cause_ranking=ranking, priority_top5=priority,
        proceed_codes=proceed_codes,
        accounting=dict(ten_groups_plus_unknown_events8=accounted_n,
            ten_groups_plus_unknown_weighted_events8=accounted_weight, exact_raw_event_accounting=True),
        design=dict(strata_n=43, masked_PSU_n=87, eligible_common_full_inventory_equal=True,
            PSU_per_stratum_distribution={'2': 42, '3': 1}),
        folds=dict(outer_salt=rules['outer_salt'], inner_salt_pattern=rules['inner_salt_pattern'],
            constructed_on_eligible_before_common=True, outer_mapping_sha256=_fingerprint(outer_map),
            outer_evaluation_PSU_counts={str(i): sum(value == i for value in outer_map.values()) for i in range(5)},
            outer=fold_diagnostics),
        weight_rule='(2/3)*WTSSCB4Y for1999/2001; (1/3)*WTSSCB2Y for2003',
        fitted_models=False, predictive_metrics_computed=False, participant_records_exported=False,
        threshold_scope=rules['threshold_scope'])
