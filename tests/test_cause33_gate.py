"""Synthetic-only cause/count qualification; no real sources, fits or metrics."""
import ast
import builtins
import copy
import hashlib
import json
from pathlib import Path
import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd

from research.mortality.cause33 import gate as g


ROOT = Path(__file__).resolve().parents[1] / 'research/mortality'
DEFINITIONS = json.loads((ROOT / 'cause33/CAUSE_DEFINITIONS.json').read_text())
FIXED_RULES = json.loads((ROOT / 'cause33/SUPPORT_RULES.json').read_text())


def rules():
    # Only the pre-existing empirical consistency counts are inapplicable to
    # clearly synthetic fixtures. Actual support minima/salts remain unchanged.
    result = copy.deepcopy(FIXED_RULES)
    result.pop('known_prior_common_n', None)
    result.pop('known_prior_common_events8', None)
    return result


def design():
    return {str(s): [1, 2, 3] if s == 1 else [1, 2] for s in range(1, 44)}


def synthetic():
    """87 public masked units; one alive plus one event for causes1..6/unit."""
    rows = []
    for st, psus in design().items():
        st = int(st); year = 1999 if st <= 13 else 2001 if st <= 28 else 2003
        for psu in psus:
            for code in range(7):
                row = dict(SEQN=len(rows) + 1, cycle=year,
                    SDDSRVYR=1 if year == 1999 else 2 if year == 2001 else 3,
                    RIDAGEYR=50, RIAGENDR=1, RIDSTATR=2, elig=1,
                    dead=int(code != 0), cause=float(code) if code else np.nan,
                    t_exam=48. if code else 177., SDMVSTRA=st, SDMVPSU=psu,
                    stratum=f'NHANES:{st}', cluster=f'{st}:{psu}',
                    WTSSCB4Y=1.5 if year != 2003 else np.nan,
                    WTSSCB2Y=3., weight=1., weight6y=1.,
                    SSCYST=1., SSB2M=2., SSCYSTL=0., SSB2ML=0.)
                row.update({name: 0. for name in g.CORE});row['male'] = 1.
                rows.append(row)
    frame = pd.DataFrame(rows)
    return frame, np.ones(len(frame), dtype=bool), dict(full_DEMO_masked_design=design(), eligible_masked_design=design())


def run(frame=None, common=None, source=None, selected_rules=None, definitions=None):
    if frame is None:
        frame, common, source = synthetic()
    return g.qualify(frame, common, source, rules() if selected_rules is None else selected_rules,
                     DEFINITIONS if definitions is None else definitions)


def reweight(frame, index, value):
    frame.loc[index, 'WTSSCB2Y'] = 3 * value
    if frame.loc[index, 'cycle'] != 2003:frame.loc[index, 'WTSSCB4Y'] = 1.5 * value
    frame.loc[index, ['weight', 'weight6y']] = value


class AccountingTests(unittest.TestCase):
    def test_ten_disjoint_groups_plus_unknown_and_common_denominator(self):
        result = run();n = 609
        self.assertEqual(result['common']['all_cause_events8'], 522)
        self.assertEqual(result['common']['all_time_deaths_n'], 522)
        self.assertEqual(result['eligible']['all_time_cause_counts']['unknown'], 0)
        self.assertEqual(result['common']['alive_n'], 87)
        self.assertEqual(result['accounting']['ten_groups_plus_unknown_events8'], 522)
        self.assertEqual(set(result['causes']), set(g.CODES))
        for code, cause in result['causes'].items():
            row = cause['overall'];events = 87 if int(code) <= 6 else 0
            self.assertEqual(row['n'], n);self.assertEqual(row['events8'], events)
            self.assertEqual(row['noncase_rows'], n - events)
            self.assertEqual(row['label_unknown_rows'], 0)
            self.assertEqual(row['common_weight_sum'], n)
            self.assertEqual(row['weighted_incidence8'], events / n)
            self.assertEqual(sum(c['events8'] for c in cause['per_cycle'].values()), events)
        self.assertEqual(result['defined_cause_ranking'], [f'{i:03d}' for i in range(1, 10)])
        self.assertEqual(result['priority_top5'], ['001', '002', '003', '004', '005'])
        self.assertEqual(result['proceed_codes'], result['priority_top5'])
        self.assertTrue(result['gate']['later_comparisons_allowed'])

    def test_competing_deaths_are_noncase_not_censored(self):
        result = run();row = result['causes']['001']['overall']
        # 87 survivors AND435 competing-cause deaths are observed noncases.
        self.assertEqual(row['noncase_rows'], 87 + 5 * 87)
        self.assertEqual(row['events8'] + row['noncase_rows'], row['n'])

    def test_zero_and_exact96_included_after96_excluded(self):
        d, m, a = synthetic();idx = d.index[d.cause.eq(1)][:3]
        d.loc[idx, 't_exam'] = [0., 96., 97.]
        result = run(d, m, a)
        self.assertEqual(result['causes']['001']['overall']['events8'], 86)
        self.assertEqual(result['common']['all_time_cause_counts']['001'], 87)
        self.assertEqual(result['causes']['001']['overall']['noncase_rows'], 609 - 86)

    def test_unknown_early_blocks_all_without_negative_unknown_labels(self):
        for time in (0., 96.):
            d, m, a = synthetic();d.loc[0, ['dead', 't_exam']] = [1., time]
            result = run(d, m, a)
            self.assertFalse(result['gate']['source_endpoint_qualified'])
            self.assertFalse(result['gate']['later_comparisons_allowed'])
            self.assertEqual(result['proceed_codes'], [])
            self.assertEqual(result['unknown_cause']['events8'], 1)
            self.assertEqual(result['unknown_cause']['event_Kish_n'], 1.)
            self.assertEqual(result['unknown_cause']['maximum_event_weight_fraction'], 1.)
            self.assertIsNone(result['unknown_cause']['noncase_rows'])
            self.assertEqual(result['accounting']['ten_groups_plus_unknown_events8'], 523)
            for group in result['causes'].values():
                row = group['overall']
                self.assertEqual(row['label_unknown_rows'], 1)
                self.assertEqual(row['noncase_rows'], row['n'] - row['events8'] - 1)
                self.assertFalse(group['proceed_to_later_protocol'])
                for cell in group['fold_cells']:
                    self.assertEqual(cell['n'], cell['events8'] + cell['noncase_rows'] + cell['label_unknown_rows'])

    def test_unknown_after96_is_known_noncase_and_reported_alltime(self):
        d, m, a = synthetic();d.loc[0, ['dead', 't_exam']] = [1., 97.]
        result = run(d, m, a)
        self.assertTrue(result['gate']['source_endpoint_qualified'])
        self.assertEqual(result['common']['unknown_cause_deaths8'], 0)
        self.assertEqual(result['common']['unknown_cause_deaths_after96'], 1)
        self.assertEqual(result['common']['all_time_cause_counts']['unknown'], 1)
        self.assertEqual(result['causes']['001']['overall']['noncase_rows'], 522)

    def test_unknown_early_only_outside_common_does_not_block(self):
        d, m, a = synthetic();d.loc[0, ['dead', 't_exam', 'sbp']] = [1., 0., np.nan];m[0] = False
        result = run(d, m, a)
        self.assertTrue(result['gate']['source_endpoint_qualified'])
        self.assertEqual(result['eligible']['unknown_cause_deaths8'], 1)
        self.assertEqual(result['common']['unknown_cause_deaths8'], 0)
        self.assertEqual(result['common']['n'], 608)

    def test_zero_and_single_event_concentration_and_weight_accounting(self):
        d, m, a = synthetic();d.loc[0, ['dead', 'cause', 't_exam']] = [1., 9., 96.];reweight(d, 0, 7.)
        result = run(d, m, a);single = result['causes']['009']['overall'];zero = result['causes']['010']['overall']
        self.assertEqual(single['weighted_events8'], 7.)
        self.assertEqual(single['event_weight_squared_sum'], 49.)
        self.assertEqual(single['event_Kish_n'], 1.)
        self.assertEqual(single['event_PSU_n'], 1);self.assertEqual(single['event_strata_n'], 1)
        self.assertEqual(single['max_event_weight8'], 7.)
        self.assertEqual(single['maximum_event_weight_fraction'], 1.)
        for field in ('event_Kish_n', 'max_event_weight8', 'maximum_event_weight_fraction'):
            self.assertIsNone(zero[field])
        self.assertEqual(zero['event_weight_squared_sum'], 0.)
        self.assertEqual(result['accounting']['ten_groups_plus_unknown_weighted_events8'], 529.)
        self.assertEqual(single['max_event_weight8'], max(c['max_event_weight8'] or 0 for c in result['causes']['009']['per_cycle'].values()))


class RankSupportTests(unittest.TestCase):
    def test_83_threshold_checks_and40_cells_with_unchanged_rules(self):
        result = run()
        for code, group in result['causes'].items():
            self.assertEqual(group['support']['threshold_check_count'], 83)
            self.assertEqual(len(group['fold_cells']), 40)
            self.assertEqual(group['support']['thresholds_passed'], int(code) <= 6)
        rare = result['causes']['009']['support']['failures']
        self.assertEqual(len(rare), 42)  # event/PSU globals +40event cells
        self.assertTrue(any(f['rule'] == 'outer_evaluation_event_rows_min' and f['required'] == 5 for f in rare))
        self.assertTrue(any(f['rule'] == 'inner_training_event_rows_min' and f['required'] == 10 for f in rare))

    def test_rank_before_support_no_backfill_and_other_never_ranked(self):
        d, m, a = synthetic();idx = d.index[d.cause.eq(1)]
        d.loc[idx[1:], ['dead', 'cause', 't_exam']] = [0., np.nan, 177.]
        reweight(d, idx[0], 1000.)
        # Other010 exceeds even the rare001 burden; it stays outside ranking.
        d.loc[0, ['dead', 'cause', 't_exam']] = [1., 10., 48.];reweight(d, 0, 2000.)
        result = run(d, m, a)
        self.assertEqual(result['priority_top5'], ['001', '002', '003', '004', '005'])
        self.assertFalse(result['causes']['001']['support']['thresholds_passed'])
        self.assertTrue(result['causes']['006']['support']['thresholds_passed'])
        self.assertEqual(result['causes']['006']['defined_cause_rank'], 6)
        self.assertEqual(result['proceed_codes'], ['002', '003', '004', '005'])
        self.assertIsNone(result['causes']['010']['defined_cause_rank'])
        self.assertNotIn('010', result['defined_cause_ranking'])
        self.assertFalse(result['causes']['010']['priority_top5'])

    def test_weight_then_raw_then_code_tie_break(self):
        d, m, a = synthetic();idx = d.index[d.cause.eq(2)][:2]
        d.loc[idx[0], ['dead', 'cause', 't_exam']] = [0., np.nan, 177.]
        reweight(d, idx[1], 2.)  # same weighted87 burden but raw86 vs001raw87
        result = run(d, m, a)
        self.assertLess(result['defined_cause_ranking'].index('001'), result['defined_cause_ranking'].index('002'))
        self.assertLess(result['defined_cause_ranking'].index('003'), result['defined_cause_ranking'].index('002'))

    def test_all_support_fail_does_not_allow_later_comparisons(self):
        r = rules();r['event_rows_min'] = 10000
        result = run(selected_rules=r)
        self.assertTrue(result['gate']['source_endpoint_qualified'])
        self.assertFalse(result['gate']['later_comparisons_allowed'])
        self.assertEqual(result['gate']['support_blockers'], ['no_supported_priority_cause'])
        self.assertEqual(result['proceed_codes'], [])


class StructureAndFoldTests(unittest.TestCase):
    def test_fractional_undeclared_or_nonfinite_causes_and_alive_cause_rejected(self):
        for value in (0., -1., 1.5, 11., np.inf):
            d, m, a = synthetic();d.loc[1, 'cause'] = value
            with self.subTest(value=value), self.assertRaises(ValueError):run(d, m, a)
        d, m, a = synthetic();d.loc[0, 'cause'] = 1.
        with self.assertRaises(ValueError):run(d, m, a)

    def test_invalid_status_cycle_time_and_types_rejected(self):
        for column, value in [('elig', 3.), ('dead', 2.), ('t_exam', -1.), ('t_exam', np.inf),
                              ('cycle', 2015.), ('SDDSRVYR', 2.), ('RIDSTATR', 1.),
                              ('RIDAGEYR', 75.5), ('SDMVPSU', 1.5), ('SEQN', 1.5)]:
            d, m, a = synthetic();d[column] = d[column].astype(float);d.loc[0, column] = value
            with self.subTest(column=column, value=value), self.assertRaises(ValueError):run(d, m, a)
        d, m, a = synthetic();d.loc[0, 't_exam'] = 95.
        with self.assertRaisesRegex(ValueError, 'Survivor'):run(d, m, a)
        for column in ('dead', 'weight', 'cause'):
            d, m, a = synthetic();d[column] = d[column].astype(object);d.loc[0, column] = '1'
            with self.subTest(column=column), self.assertRaisesRegex(ValueError, 'Nonnumeric'):run(d, m, a)
        d, m, a = synthetic();d['dead'] = d.dead.astype(bool)
        with self.assertRaises(ValueError):run(d, m, a)

    def test_exact_weights_and_expected_late_W4_missingness(self):
        for column, value in [('weight', 1.0000000000001), ('weight', 0.), ('weight', -1.),
                              ('WTSSCB4Y', np.nan), ('WTSSCB2Y', np.inf), ('weight6y', 2.)]:
            d, m, a = synthetic();d.loc[0, column] = value
            with self.subTest(column=column), self.assertRaises(ValueError):run(d, m, a)
        d, m, a = synthetic();late = d.index[d.cycle.eq(2003)][0];d.loc[late, 'WTSSCB4Y'] = 1.
        with self.assertRaises(ValueError):run(d, m, a)

    def test_mask_is_exact_boolean_outcome_free_domain_and_preserves_allPSUs(self):
        d, m, a = synthetic();m[0] = False
        with self.assertRaisesRegex(ValueError, 'Common mask changed'):run(d, m, a)
        with self.assertRaises(ValueError):run(d, m.astype(int), a)
        d.loc[0, 'sbp'] = np.nan
        self.assertEqual(run(d, m, a)['common']['n'], 608)
        cluster = d.cluster.eq('1:3');d.loc[cluster, 'sbp'] = np.nan;m = np.isfinite(d[list(g.CORE)]).all(axis=1).to_numpy()
        with self.assertRaisesRegex(ValueError, 'inventory'):run(d, m, a)

    def test_duplicate_fields_ids_design_labels_and_audit_rejected(self):
        d, m, a = synthetic();d.loc[0, 'SEQN'] = d.loc[1, 'SEQN']
        with self.assertRaises(ValueError):run(d, m, a)
        d, m, a = synthetic();d = pd.concat([d, d[['cause']]], axis=1)
        with self.assertRaises(ValueError):run(d, m, a)
        for column, value in [('cluster', '1'), ('stratum', '1'), ('SDMVSTRA', 29.)]:
            d, m, a = synthetic();d.loc[0, column] = value
            with self.subTest(column=column), self.assertRaises(ValueError):run(d, m, a)
        d, m, a = synthetic();a['full_DEMO_masked_design']['1'] = [1, 2]
        with self.assertRaises(ValueError):run(d, m, a)

    def test_optional_prior_count_gates_and_no_changed_cause_registry(self):
        r = rules();r['known_prior_common_n'] = 6044
        with self.assertRaisesRegex(ValueError, 'Prior common'):run(selected_rules=r)
        r['known_prior_common_n'] = 609;r['known_prior_common_events8'] = 522
        self.assertEqual(run(selected_rules=r)['common']['n'], 609)
        defs = copy.deepcopy(DEFINITIONS);defs['supported_baseline_start_years'] = [2015, 2017]
        with self.assertRaises(ValueError):run(definitions=defs)
        defs = copy.deepcopy(DEFINITIONS);defs['groups']['001']['id'] = 'all_CVD'
        with self.assertRaises(ValueError):run(definitions=defs)
        defs = copy.deepcopy(DEFINITIONS);defs['groups'] = None
        with self.assertRaises(ValueError):run(definitions=defs)
        r = rules();r['event_rows_min'] = True
        with self.assertRaises(ValueError):run(selected_rules=r)

    def test_outcome_independent_wholePSU_maps_equal_retained31_pure_functions(self):
        d, m, a = synthetic();result = run(d, m, a)
        keys = d.cluster.to_numpy();outer = g._partition_map(keys, 'b2m31:outer', 5)
        # Compile ONLY the two retained pure partition functions, importing no
        # model/analyze module and executing no numerical training primitive.
        def pure(path, name):
            tree = ast.parse(path.read_text());node = next(x for x in tree.body if isinstance(x, ast.FunctionDef) and x.name == name)
            unit = ast.Module(body=[node], type_ignores=[]);namespace = {'np': np, 'hashlib': hashlib}
            exec(compile(unit, '<retained-pure-partition>', 'exec'), namespace)
            return namespace[name]
        old_outer = pure(ROOT / 'b2m31/analyze.py', 'assign_outer')
        old_inner = pure(ROOT / 'binary28/model.py', 'psu_partition')
        np.testing.assert_array_equal(old_outer(keys), np.array([outer[k] for k in keys]))
        for fold in range(5):
            train = np.array([outer[k] != fold for k in keys]);inner = g._partition_map(keys[train], f'b2m31:inner:{fold}', 3)
            np.testing.assert_array_equal(old_inner(keys[train], f'b2m31:inner:{fold}', folds=3), [inner[k] for k in keys[train]])
        changed = d.iloc[::-1].copy();changed['dead'] = 0.;changed['cause'] = np.nan;changed['t_exam'] = 177.
        second = run(changed, np.ones(len(changed), dtype=bool), a)
        self.assertEqual(result['folds'], second['folds'])
        self.assertEqual(result['design']['masked_PSU_n'], 87)
        self.assertEqual(result['design']['PSU_per_stratum_distribution'], {'2': 42, '3': 1})

    def test_no_fitting_imports_runtime_and_no_participant_output_or_mutation(self):
        tree = ast.parse(Path(g.__file__).read_text())
        imported = {node.module.split('.')[0] for node in ast.walk(tree) if isinstance(node, ast.ImportFrom)}
        imported |= {alias.name.split('.')[0] for node in ast.walk(tree) if isinstance(node, ast.Import) for alias in node.names}
        self.assertLessEqual(imported, {'__future__', 'hashlib', 'json', 'numpy', 'pandas'})
        real_import = builtins.__import__
        def guarded(name, *args, **kwargs):
            if any(part in {'sklearn', 'scipy', 'model', 'metrics', 'analyze'} for part in name.split('.')):
                raise AssertionError('Fitting/metrics import in count gate')
            return real_import(name, *args, **kwargs)
        d, m, a = synthetic();before = d.copy(deep=True);old_mask = m.copy();old_audit = copy.deepcopy(a)
        with patch('builtins.__import__', side_effect=guarded):result = run(d, m, a)
        pd.testing.assert_frame_equal(d, before);np.testing.assert_array_equal(m, old_mask);self.assertEqual(a, old_audit)
        text = json.dumps(result, allow_nan=False)
        self.assertNotIn('SEQN', text);self.assertNotIn('predictions', text)
        self.assertFalse(result['fitted_models']);self.assertFalse(result['predictive_metrics_computed'])
        self.assertFalse(result['participant_records_exported'])


if __name__ == '__main__':
    unittest.main()
