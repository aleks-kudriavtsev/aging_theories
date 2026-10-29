"""Synthetic-only33B source/domain/endpoint tests; no real data or fitting."""
import copy
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd

from research.mortality.cause33b import cohort as c


PRIOR = Path(__file__).resolve().parents[1] / 'research/mortality/cause33'
DEFINITIONS = json.loads((PRIOR / 'CAUSE_DEFINITIONS.json').read_text())
RULES = json.loads((PRIOR / 'SUPPORT_RULES.json').read_text())


def synthetic():
    """All87 public masked units with two supported, ordered cause targets."""
    rows = []
    design = {str(s): [1, 2, 3] if s == 1 else [1, 2] for s in range(1, 44)}
    for stratum, psus in design.items():
        st = int(stratum); year = 1999 if st <= 13 else 2001 if st <= 28 else 2003
        for psu in psus:
            for code in (0, 1, 2):
                weight = 2. if code == 2 else 1.
                row = dict(SEQN=99888000 + len(rows), cycle=year,
                    SDDSRVYR=1 if year == 1999 else 2 if year == 2001 else 3,
                    RIDAGEYR=50, RIAGENDR=1, RIDSTATR=2, elig=1,
                    dead=int(code != 0), cause=float(code) if code else np.nan,
                    t_exam=48. if code else 177., SDMVSTRA=st, SDMVPSU=psu,
                    stratum=f'NHANES:{st}', cluster=f'{st}:{psu}',
                    WTSSCB4Y=1.5 * weight if year != 2003 else np.nan,
                    WTSSCB2Y=3. * weight, weight=weight, weight6y=weight,
                    SSCYST=1., SSB2M=2., SSCYSTL=0., SSB2ML=0.)
                row.update({name: 0. for name in c.gate.CORE}); row['male'] = 1.
                rows.append(row)
    rules = copy.deepcopy(RULES)
    rules.pop('known_prior_common_n'); rules.pop('known_prior_common_events8')
    return (pd.DataFrame(rows), np.ones(len(rows), dtype=bool),
            dict(full_DEMO_masked_design=design, eligible_masked_design=design), rules)


def frozen_result(result):
    frozen = copy.deepcopy(result)
    frozen.update(schema='cause33-qualified-counts/1', plan_commit='a' * 40,
        lock_sha256='b' * 64, runtime_python='3.synthetic',
        prior_cohort_use='all_cycles_previously_used_no_independent_validation',
        comparison_models_fitted=False, independent_binary_probabilities_summed=False,
        joint_integration_and_Fine_Gray='deferred_not_implemented_or_run',
        historical_scientific_files_verified=99, scientific_python_sha256={})
    return frozen


class ProjectionTests(unittest.TestCase):
    def setUp(self):
        self.frame, self.common, self.audit, self.rules = synthetic()
        self.result = c.gate.qualify(self.frame, self.common, self.audit, self.rules, DEFINITIONS)
        self.frozen = frozen_result(self.result)

    def qualify(self, frozen=None):
        return c.reproduce_gate(self.frame, self.common, self.audit, self.rules,
            DEFINITIONS, self.frozen if frozen is None else frozen)

    def test_entire_scientific_gate_reproduces_with_counts_weights_folds_support(self):
        before = self.frame.copy(deep=True); mask = self.common.copy()
        result, proof = self.qualify()
        self.assertEqual(result, self.result)
        self.assertEqual(result['proceed_codes'], ['002', '001'])
        self.assertTrue(proof['exact_structure_and_counts'])
        self.assertEqual(proof['numeric_tolerance'], dict(rtol=1e-12, atol=1e-9))
        self.assertGreater(proof['floating_values_checked'], 100)
        self.assertEqual(proof['floating_values_different'], 0)
        self.assertEqual(proof['recomputed_scientific_sha256'], proof['frozen_scientific_sha256'])
        pd.testing.assert_frame_equal(self.frame, before)
        np.testing.assert_array_equal(self.common, mask)

    def test_known_runner_provenance_may_change_without_dropping_scientific_fields(self):
        self.frozen.update(runtime_python='different platform', scientific_python_sha256={'new.py': 'c' * 64})
        _, proof = self.qualify()
        self.assertEqual(proof['floating_values_different'], 0)
        for change in ('unknown-extra', 'missing-scientific'):
            frozen = copy.deepcopy(self.frozen)
            if change == 'unknown-extra': frozen['new_unreviewed_field'] = 1
            else: frozen.pop('common')
            with self.subTest(change=change), self.assertRaisesRegex(ValueError, 'field set'):
                self.qualify(frozen)

    def test_counts_fold_fingerprints_ranking_support_and_weight_changes_block(self):
        changes = [
            lambda f: f['common'].__setitem__('n', f['common']['n'] + 1),
            lambda f: f['folds'].__setitem__('outer_mapping_sha256', 'f' * 64),
            lambda f: f['priority_top5'].reverse(),
            lambda f: f['causes']['001']['support'].__setitem__('thresholds_passed', False),
            lambda f: f['per_cycle']['1999']['common'].__setitem__('weight_sum', 999.),
        ]
        for index, change in enumerate(changes):
            frozen = copy.deepcopy(self.frozen); change(frozen)
            with self.subTest(change=index), self.assertRaisesRegex(ValueError, '33A scientific'):
                self.qualify(frozen)

    def test_only_finite_float_roundoff_has_the_locked_tolerance(self):
        path = self.frozen['common']
        path['weight_sum'] = float(np.nextafter(path['weight_sum'], np.inf))
        _, proof = self.qualify()
        self.assertEqual(proof['floating_values_different'], 1)
        self.assertGreater(proof['maximum_float_abs_difference'], 0.)
        self.assertNotEqual(proof['recomputed_scientific_sha256'], proof['frozen_scientific_sha256'])
        for value in (np.inf, np.nan, self.result['common']['weight_sum'] + 1e-5):
            frozen = copy.deepcopy(self.frozen); frozen['common']['weight_sum'] = value
            with self.subTest(value=value), self.assertRaisesRegex(ValueError, 'floating aggregate'):
                self.qualify(frozen)

    def test_integer_counts_and_booleans_cannot_be_floats_or_each_other(self):
        for key, value in [('n', float(self.result['common']['n'])), ('n', True)]:
            frozen = copy.deepcopy(self.frozen); frozen['common'][key] = value
            with self.subTest(value=value), self.assertRaisesRegex(ValueError, 'type changed'):
                self.qualify(frozen)
        frozen = copy.deepcopy(self.frozen); frozen['fitted_models'] = 0
        with self.assertRaisesRegex(ValueError, 'type changed'): self.qualify(frozen)

    def test_unknown_early_cause_stops_before_any_reproduction_or_negative_label(self):
        self.frame.loc[0, ['dead', 't_exam']] = [1., 96.]
        with self.assertRaisesRegex(ValueError, 'blocks all comparisons'):
            self.qualify()

    def test_supported_causes_are_fixed_not_any_supported_subset_or_reordering(self):
        for codes in (['001', '002'], ['002'], ['002', '001', '003']):
            result = copy.deepcopy(self.result); result['proceed_codes'] = codes
            with self.subTest(codes=codes), patch.object(c.gate, 'qualify', return_value=result), \
                    self.assertRaisesRegex(ValueError, 'supported target'):
                self.qualify()


class EndpointTests(unittest.TestCase):
    def frame(self):
        return pd.DataFrame(dict(dead=[1., 1., 1., 1., 1., 1., 0., 0., 1., 1.],
            cause=[2., 1., 3., 10., 2., 1., np.nan, np.nan, np.nan, np.nan],
            t_exam=[0., 96., 40., 96., 97., 150., 96., 177., 110., 12.]),
            index=range(100, 110))

    def test_common_order_exact_horizon_competing_deaths_and_late_unknown(self):
        frame = self.frame(); common = np.ones(len(frame), dtype=bool); common[-1] = False
        before = frame.copy(deep=True)
        labels = c.cause_targets(frame, common)
        self.assertEqual(list(labels), ['002', '001'])
        for code, index in [('002', 0), ('001', 1)]:
            expected = np.zeros(9, dtype=bool); expected[index] = True
            np.testing.assert_array_equal(labels[code], expected)
            self.assertEqual(labels[code].dtype, bool)
        labels['002'][0] = False
        pd.testing.assert_frame_equal(frame, before)
        self.assertTrue(c.cause_targets(frame, common)['002'][0])

    def test_unknown_common_at_zero_and96_stops_all_targets_without_exclusion(self):
        for month in (0., 96.):
            frame = self.frame(); frame.loc[109, 't_exam'] = month
            with self.subTest(month=month), self.assertRaisesRegex(ValueError, 'Unknown common early'):
                c.cause_targets(frame, np.ones(len(frame), dtype=bool))

    def test_invalid_status_cause_and_survivor_followup_are_not_repaired(self):
        changes = [('dead', .5), ('dead', np.nan), ('cause', .5), ('cause', 11.),
                   ('cause', np.inf), ('t_exam', -1.), ('t_exam', np.nan)]
        for field, value in changes:
            frame = self.frame(); frame.loc[100, field] = value
            with self.subTest(field=field, value=value), self.assertRaisesRegex(ValueError, 'Invalid33B'):
                c.cause_targets(frame, np.ones(len(frame), dtype=bool))
        for field, value in [('cause', 1.), ('t_exam', 95.)]:
            frame = self.frame(); frame.loc[106, field] = value
            with self.subTest(field=field), self.assertRaisesRegex(ValueError, 'Invalid33B'):
                c.cause_targets(frame, np.ones(len(frame), dtype=bool))

    def test_boolean_mask_index_and_numeric_endpoint_types_are_strict(self):
        frame = self.frame(); common = np.ones(len(frame), dtype=bool); common[-1] = False
        for wrong in (common.astype(int), common[:-1], pd.Series(common),
                      pd.Series(common, index=frame.index[::-1])):
            with self.subTest(mask=type(wrong).__name__), self.assertRaisesRegex(ValueError, 'boolean common'):
                c.cause_targets(frame, wrong)
        nullable = pd.Series(common, index=frame.index, dtype='boolean')
        nullable.iloc[0] = pd.NA
        with self.assertRaisesRegex(ValueError, 'boolean common'):
            c.cause_targets(frame, nullable)
        for dtype in (str, bool):
            changed = frame.copy(); changed['dead'] = changed.dead.astype(dtype)
            with self.subTest(dtype=dtype), self.assertRaisesRegex(ValueError, 'Numeric endpoint'):
                c.cause_targets(changed, common)


class SourceAdapterTests(unittest.TestCase):
    def setUp(self):
        self.frame, self.common, self.inherited, self.rules = synthetic()
        result = c.gate.qualify(self.frame, self.common, self.inherited, self.rules, DEFINITIONS)
        self.frozen = frozen_result(result)
        self.base = dict(schema='synthetic-manifest/1', files=[])
        self.wrapper = dict(schema='synthetic33b-wrapper/1', source='synthetic-only')
        self.directory = tempfile.TemporaryDirectory(); self.addCleanup(self.directory.cleanup)
        self.root = Path(self.directory.name); self.prior = self.root / 'prior'; self.prior.mkdir()
        (self.root / 'CONTRACT.json').write_text(json.dumps(dict(
            source_projection_float_rtol=1e-12, source_projection_float_atol=1e-9)))
        for name, value in [('GATE_RESULTS.json', self.frozen), ('SUPPORT_RULES.json', self.rules),
                            ('CAUSE_DEFINITIONS.json', DEFINITIONS)]:
            (self.prior / name).write_text(json.dumps(value))
        for name, value in [('ROOT', self.root), ('PRIOR', self.prior)]:
            patcher = patch.object(c, name, value); patcher.start(); self.addCleanup(patcher.stop)
        loader = patch.object(c, '_load_prior', return_value=(self.frame, self.common,
            self.common.copy(), self.inherited, self.base))
        self.loader = loader.start(); self.addCleanup(loader.stop)
        verify = patch.object(c, '_verify_wrapper', return_value=(self.base, self.wrapper))
        self.verify = verify.start(); self.addCleanup(verify.stop)

    def test_source_verified_before_parsing_and_rechecked_before_aggregate_return(self):
        order = []
        self.verify.side_effect = lambda source: (order.append('wrapper') or (self.base, self.wrapper))
        self.loader.side_effect = lambda source: (order.append('parse') or
            (self.frame, self.common, self.common.copy(), self.inherited, self.base))
        frame, common, audit, manifest = c.load(self.root)
        self.assertEqual(order, ['wrapper', 'parse', 'wrapper'])
        self.assertIs(frame, self.frame); self.assertIs(common, self.common)
        self.assertEqual(audit['reproduced33A_gate']['proceed_codes'], ['002', '001'])
        self.assertTrue(audit['scientific_gate_equality']['exact_structure_and_counts'])
        self.assertEqual(manifest['source'], self.base); self.assertEqual(manifest['intake'], self.wrapper)
        self.assertEqual(audit['frozen33A_gate_sha256'], manifest['frozen33A_gate_sha256'])
        self.assertFalse(audit['participant_records_exported']); self.assertFalse(audit['labels_exported'])
        serialized = json.dumps(dict(audit=audit, manifest=manifest), allow_nan=False)
        self.assertNotIn('SEQN', serialized); self.assertNotIn('99888000', serialized)

    def test_wrapper_failure_prevents_any_source_parsing(self):
        self.verify.side_effect = ValueError('fixed wrapper blocked')
        with self.assertRaisesRegex(ValueError, 'fixed wrapper blocked'): c.load(self.root)
        self.loader.assert_not_called()

    def test_contract_tolerance_changes_cannot_widen_the_scientific_gate(self):
        contract = self.root / 'CONTRACT.json'
        contract.write_text(json.dumps(dict(source_projection_float_rtol=1e-3,
                                           source_projection_float_atol=1e-9)))
        with self.assertRaisesRegex(ValueError, 'tolerance differs'): c.load(self.root)
        self.loader.assert_not_called()

    def test_changed_base_manifest_stops_before_gate(self):
        self.loader.return_value = (self.frame, self.common, None, self.inherited, {'changed': True})
        with patch.object(c.gate, 'qualify', side_effect=AssertionError('must not qualify')):
            with self.assertRaisesRegex(ValueError, 'manifest changed'): c.load(self.root)

    def test_changed_wrapper_during_parsing_blocks_return(self):
        self.verify.side_effect = [(self.base, self.wrapper), (self.base, {'changed': True})]
        with self.assertRaisesRegex(ValueError, 'provenance changed'): c.load(self.root)

    def test_changed_frozen_reference_during_qualification_blocks_return(self):
        reproduce = c.reproduce_gate
        def changed(*args):
            result = reproduce(*args)
            with (self.prior / 'GATE_RESULTS.json').open('a') as stream: stream.write('\n')
            return result
        with patch.object(c, 'reproduce_gate', side_effect=changed):
            with self.assertRaisesRegex(ValueError, 'provenance changed'): c.load(self.root)


if __name__ == '__main__':
    unittest.main()
