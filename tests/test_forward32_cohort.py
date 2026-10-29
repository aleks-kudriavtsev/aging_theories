"""Synthetic temporal/source gates; no real participant input or model fits."""
import copy
import hashlib
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd

from research.mortality.forward32 import cohort as c
from research.mortality.forward32 import retrieve as r


LOCK = dict(minimum_positive_weight_event_rows_per_part=20,
            minimum_positive_weight_nonevent_rows_per_part=20,
            inner_salt='forward32:inner:v1')


def full_design():
    return {str(st): [1, 2, 3] if st == 1 else [1, 2] for st in range(1, 44)}


def synthetic():
    """All43 public masked strata/87PSUs; one event and nonevent perPSU."""
    rows = []
    for st, psus in full_design().items():
        st = int(st)
        cycle = 1999 if st <= 13 else 2001 if st <= 28 else 2003
        for psu in psus:
            for event in (0, 1):
                seq = len(rows) + 1
                w4 = 10. + st / 10 + psu if cycle != 2003 else np.nan
                w2 = 20. + st / 5 + psu
                row = dict(SEQN=seq, cycle=cycle, SDMVSTRA=st, SDMVPSU=psu,
                    stratum=str(st), cluster=f'{st}:{psu}',
                    WTSSCB4Y=w4, WTSSCB2Y=w2,
                    weight=(2 / 3) * w4 if cycle != 2003 else (1 / 3) * w2,
                    dead=event, t_exam=48. if event else 177., SSPRIS=0.,
                    SSCYSTL=0., SSB2ML=0., SSCYST=1 + seq / 1000,
                    SSB2M=2 + seq / 1000)
                row.update({name: (j + 1) * seq / 100 for j, name in enumerate(c.CORE)})
                rows.append(row)
    frame = pd.DataFrame(rows)
    return frame, np.ones(len(frame), dtype=bool), full_design()


def qualify(frame=None, common=None, design=None):
    if frame is None:
        frame, default_common, default_design = synthetic()
        common = default_common if common is None else common
        design = default_design if design is None else design
    return c.temporal_qualification(frame, common, design, LOCK)


def wrapper_fixture(directory):
    """Only synthetic JSON metadata; no source/XPT bytes are created."""
    directory = Path(directory)
    base = dict(schema='synthetic-base-manifest/1', files=[
        dict(file='1999_SSCARD_A.xpt', sha256='a' * 64)])
    (directory / 'manifest.json').write_text(json.dumps(base))
    wrapper = dict(schema='forward32-retrieval/1', plan_commit=r.PLAN_COMMIT,
        lock_sha256=r.sha(r.ROOT / 'LOCK.json'),
        base_manifest_sha256=r.sha(directory / 'manifest.json'),
        source_fingerprints={'1999_SSCARD_A.xpt': 'a' * 64})
    (directory / 'forward32_manifest.json').write_text(json.dumps(wrapper))
    return base, wrapper


class NativeWeightTests(unittest.TestCase):
    def test_native_assignment_is_w4_early_w2_late_not_original_w6(self):
        frame, _, _ = synthetic(); before = frame.copy(deep=True)
        d = c.native_period_frame(frame)
        early = frame.cycle.isin([1999, 2001]); late = ~early
        np.testing.assert_array_equal(d.loc[early, 'weight'], frame.loc[early, 'WTSSCB4Y'])
        np.testing.assert_array_equal(d.loc[late, 'weight'], frame.loc[late, 'WTSSCB2Y'])
        np.testing.assert_array_equal(d.weight6y_original, frame.weight)
        np.testing.assert_allclose(d.loc[early, 'weight'], frame.loc[early, 'weight'] * 1.5)
        np.testing.assert_allclose(d.loc[late, 'weight'], frame.loc[late, 'weight'] * 3)
        self.assertEqual(d.partition.tolist(), np.where(early, 'train', 'late').tolist())
        self.assertTrue(frame.loc[late, 'WTSSCB4Y'].isna().all())
        self.assertEqual(len(d), len(frame)); pd.testing.assert_frame_equal(frame, before)
        self.assertIsNot(d, frame)

    def test_period_weight_sums_raw_sources_effective_n_and_common_coverage(self):
        frame, common, design = synthetic()
        common[0] = False; common[np.flatnonzero(frame.cycle.eq(2003))[0]] = False
        before = frame.copy(deep=True); before_common = common.copy()
        d, audit = qualify(frame, common, design)
        self.assertTrue(audit['gate']['temporal_evaluation_allowed'])
        for part, mask, source in [('train', frame.cycle.ne(2003), 'WTSSCB4Y'),
                                   ('late', frame.cycle.eq(2003), 'WTSSCB2Y')]:
            p = audit['periods'][part]; eligible = frame.loc[mask]
            target = frame.loc[mask.to_numpy() & common]
            w = target[source].to_numpy()
            self.assertAlmostEqual(p['eligible']['native_weight_sum'], eligible[source].sum())
            self.assertAlmostEqual(p['common']['native_weight_sum'], w.sum())
            self.assertAlmostEqual(p['common']['native_weight_squared_sum'], w @ w)
            self.assertAlmostEqual(p['common']['effective_n'], w.sum() ** 2 / (w @ w))
            self.assertAlmostEqual(p['common_native_weighted_coverage'], w.sum() / eligible[source].sum())
            self.assertEqual(p['common']['selected_native_source_max_abs_difference'], 0.)
            self.assertAlmostEqual(p['common']['original_W6_sum'], target.weight.sum())
            self.assertEqual(p['weight_variable'], source)
            for cycle, group in target.groupby('cycle'):
                sums = p['common']['per_cycle_weights'][str(cycle)]
                self.assertAlmostEqual(sums['native_weight_sum'], group[source].sum())
                self.assertAlmostEqual(sums['original_W6_sum'], group.weight.sum())
                self.assertAlmostEqual(sums['source_W2_sum'], group.WTSSCB2Y.sum())
        self.assertEqual(audit['periods']['late']['eligible']['source_W4_missing_n'], 60)
        self.assertEqual(audit['periods']['late']['eligible']['n'], 60)
        self.assertTrue(audit['no_six_year_pooled_score']); self.assertTrue(audit['no_late_recalibration'])
        self.assertTrue(audit['gate']['support_threshold_is_not_power_or_calibration_certificate'])
        self.assertFalse(audit['assay_batch_individually_identified'])
        pd.testing.assert_frame_equal(frame, before)
        np.testing.assert_array_equal(common, before_common)
        self.assertEqual(len(d), len(frame))

    def test_irrelevant_native_source_missingness_does_not_drop_other_period(self):
        frame, common, design = synthetic()
        frame.loc[frame.cycle.ne(2003), 'WTSSCB2Y'] = np.nan
        d, audit = qualify(frame, common, design)
        self.assertEqual(len(d), 174)
        self.assertTrue(audit['gate']['temporal_evaluation_allowed'])
        self.assertTrue(np.isfinite(d.weight).all())

    def test_selected_native_zero_negative_or_missing_blocks_and_diagnostic_is_finite_json(self):
        for part, source in [('train', 'WTSSCB4Y'), ('late', 'WTSSCB2Y')]:
            for value in (0., -1., np.nan, np.inf):
                frame, common, design = synthetic()
                mask = frame.cycle.ne(2003) if part == 'train' else frame.cycle.eq(2003)
                frame.loc[frame.index[mask][0], source] = value
                with self.subTest(part=part, value=value):
                    _, audit = qualify(frame, common, design)
                    self.assertIn(part + ':native_weight_unusable', audit['gate']['blockers'])
                    self.assertFalse(audit['gate']['temporal_evaluation_allowed'])
                    json.dumps(audit, allow_nan=False)

    def test_unknown_or_missing_cycle_is_rejected_instead_of_becoming_late(self):
        for cycle in (1998, 2000, 2002, 2005, np.nan):
            frame, _, _ = synthetic(); frame.loc[0, 'cycle'] = cycle
            with self.subTest(cycle=cycle), self.assertRaisesRegex(ValueError, 'Unsupported temporal cycle'):
                c.native_period_frame(frame)


class DesignAndSupportTests(unittest.TestCase):
    def test_canonical43_strata87_psus_and_psu3_preserved_without_false_overlap(self):
        _, audit = qualify()
        self.assertTrue(audit['gate']['temporal_evaluation_allowed'])
        self.assertEqual(audit['participant_overlap_n'], 0)
        self.assertEqual(audit['masked_PSU_overlap_n'], 0)
        self.assertEqual(audit['periods']['train']['common']['strata_n'], 28)
        self.assertEqual(audit['periods']['train']['common']['psu_n'], 57)
        self.assertEqual(audit['periods']['late']['common']['strata_n'], 15)
        self.assertEqual(audit['periods']['late']['common']['psu_n'], 30)
        self.assertEqual(audit['periods']['train']['eligible_masked_design']['1'], [1, 2, 3])
        self.assertEqual(audit['periods']['train']['common_masked_design']['1'], [1, 2, 3])
        # Reused PSU labels1/2 across strata are not shared masked clusters.
        self.assertEqual(audit['periods']['late']['eligible_masked_design']['29'], [1, 2])

    def test_lost_eligible_psu_third_psu_or_whole_stratum_blocks_without_relabeling(self):
        for st, psu in ((1, 3), (2, 2), (29, None)):
            frame, _, design = synthetic()
            remove = frame.SDMVSTRA.eq(st)
            if psu is not None: remove &= frame.SDMVPSU.eq(psu)
            frame = frame.loc[~remove].copy(); common = np.ones(len(frame), dtype=bool)
            part = 'train' if st <= 28 else 'late'
            with self.subTest(st=st, psu=psu):
                _, audit = qualify(frame, common, design)
                self.assertFalse(audit['gate']['temporal_evaluation_allowed'])
                self.assertIn(part + ':full_masked_design_inventory_lost', audit['gate']['blockers'])
                if st == 2: self.assertIn('train:singleton_stratum', audit['gate']['blockers'])
                if st == 1:
                    self.assertEqual(audit['periods']['train']['eligible_masked_design']['1'], [1, 2])
                    self.assertNotIn('train:singleton_stratum', audit['gate']['blockers'])

    def test_complete_case_loss_of_any_entire_psu_blocks_and_keeps_eligible_reference(self):
        for part, st, psu in [('train', 1, 3), ('late', 29, 2)]:
            frame, common, design = synthetic()
            common[(frame.SDMVSTRA.eq(st) & frame.SDMVPSU.eq(psu)).to_numpy()] = False
            _, audit = qualify(frame, common, design)
            with self.subTest(part=part):
                self.assertFalse(audit['gate']['temporal_evaluation_allowed'])
                self.assertIn(part + ':common_domain_lost_eligible_PSU_fold_equivalence',
                              audit['gate']['blockers'])
                self.assertEqual(audit['periods'][part]['eligible_masked_design'][str(st)], design[str(st)])
                self.assertEqual(audit['periods'][part]['common_masked_design'][str(st)],
                                 [p for p in design[str(st)] if p != psu])
                self.assertEqual(audit['periods'][part]['common']['psu_n'],
                                 audit['periods'][part]['eligible']['psu_n'] - 1)

    def test_partial_row_filter_preserves_eligible_fold_mapping_and_independent_hash(self):
        frame, common, design = synthetic()
        common[0] = False; common[np.flatnonzero(frame.cycle.eq(2003))[0]] = False
        _, audit = qualify(frame, common, design)
        self.assertTrue(audit['gate']['temporal_evaluation_allowed'])
        self.assertTrue(audit['common_eligible_PSU_fold_equivalence_required'])
        for part in ('train', 'late'):
            self.assertEqual(audit['periods'][part]['common_masked_design'],
                             audit['periods'][part]['eligible_masked_design'])
        early = frame.cycle.ne(2003).to_numpy()
        eligible_keys = sorted(set(frame.loc[early, 'cluster']))
        hash_order = sorted(eligible_keys, key=lambda key: hashlib.sha256(
            (LOCK['inner_salt'] + ':' + key).encode()).hexdigest())
        expected = {key: hash_order.index(key) % 3 for key in eligible_keys}
        exported = audit['eligible_derived_early_inner_fold_mapping']
        self.assertEqual(exported, expected)
        self.assertEqual(len(exported), 57)
        common_clusters = frame.loc[early & common, 'cluster'].to_numpy()
        common_folds = c.psu_partition(common_clusters, LOCK['inner_salt'], folds=3)
        common_mapping = {key: int(common_folds[np.flatnonzero(common_clusters == key)[0]])
                          for key in sorted(set(common_clusters))}
        self.assertEqual(common_mapping, exported)
        def mapping_hash(mapping):
            return hashlib.sha256(json.dumps(mapping, sort_keys=True, separators=(',', ':')).encode()).hexdigest()
        self.assertEqual(mapping_hash(common_mapping), mapping_hash(expected))
        self.assertEqual(mapping_hash(exported), mapping_hash(expected))
        # Outcomes, weights and row order cannot enter the eligible identity hash.
        changed = frame.iloc[::-1].copy(); changed['dead'] = 1 - changed.dead
        changed['t_exam'] = np.where(changed.dead.eq(1), 48., 177.)
        changed.loc[changed.cycle.ne(2003), 'WTSSCB4Y'] *= 7.
        _, changed_audit = qualify(changed, common[::-1], design)
        self.assertEqual(changed_audit['eligible_derived_early_inner_fold_mapping'], expected)

    def test_incomplete_or_extra_full_design_reference_fails_explicitly(self):
        for mutation in ('missing', 'extra'):
            frame, common, design = synthetic()
            if mutation == 'missing': del design['43']
            else: design['44'] = [1, 2]
            with self.subTest(mutation=mutation), self.assertRaisesRegex(ValueError, 'design reference'):
                qualify(frame, common, design)

    def test_participant_or_masked_cluster_overlap_is_blocked(self):
        for column in ('SEQN', 'cluster'):
            frame, common, design = synthetic()
            late_row = frame.index[frame.cycle.eq(2003)][0]
            frame.loc[late_row, column] = frame.loc[0, column]
            _, audit = qualify(frame, common, design)
            self.assertIn('temporal_participant_or_masked_PSU_overlap', audit['gate']['blockers'])
            self.assertGreater(audit['participant_overlap_n'] if column == 'SEQN' else audit['masked_PSU_overlap_n'], 0)

    def test_cycle_stratum_disagreement_cannot_silently_move_inventory(self):
        frame, common, design = synthetic()
        frame.loc[frame.SDMVSTRA.eq(29), 'cycle'] = 2001
        frame.loc[frame.SDMVSTRA.eq(29), 'WTSSCB4Y'] = 10.
        _, audit = qualify(frame, common, design)
        self.assertIn('train:full_masked_design_inventory_lost', audit['gate']['blockers'])
        self.assertIn('late:full_masked_design_inventory_lost', audit['gate']['blockers'])

    def test_both_classes_support_boundary20_passes19_blocks_for_each_period(self):
        for part in ('train', 'late'):
            for support_class in (0, 1):
                for count in (20, 19):
                    frame, common, design = synthetic()
                    mask = frame.cycle.ne(2003) if part == 'train' else frame.cycle.eq(2003)
                    indices = frame.index[mask]
                    frame.loc[indices, 'dead'] = 1 - support_class
                    frame.loc[indices[:count], 'dead'] = support_class
                    frame['t_exam'] = np.where(frame.dead.eq(1), 48., 177.)
                    _, audit = qualify(frame, common, design)
                    key = 'event' if support_class else 'nonevent'
                    blocker = part + ':' + key + '_support_below_predeclared_minimum'
                    with self.subTest(part=part, support_class=support_class, count=count):
                        self.assertEqual(audit['periods'][part]['common']['positive_weight_' + key + 's8'], count)
                        self.assertEqual(blocker in audit['gate']['blockers'], count == 19)
                        self.assertEqual(audit['gate']['temporal_evaluation_allowed'], count == 20)

    def test_zero_weight_event_does_not_count_toward_twenty(self):
        frame, common, design = synthetic(); late = frame.index[frame.cycle.eq(2003)]
        frame.loc[late, 'dead'] = 0; frame.loc[late[:20], 'dead'] = 1
        frame['t_exam'] = np.where(frame.dead.eq(1), 48., 177.)
        frame.loc[late[0], 'WTSSCB2Y'] = 0.
        _, audit = qualify(frame, common, design)
        self.assertEqual(audit['periods']['late']['common']['events8'], 20)
        self.assertEqual(audit['periods']['late']['common']['positive_weight_events8'], 19)
        self.assertIn('late:event_support_below_predeclared_minimum', audit['gate']['blockers'])

    def test_early_censor_gate_includes_eligible_rows_outside_common_domain(self):
        for part in ('train', 'late'):
            frame, common, design = synthetic()
            mask = frame.cycle.ne(2003) if part == 'train' else frame.cycle.eq(2003)
            i = frame.index[mask & frame.dead.eq(0)][0]
            frame.loc[i, 't_exam'] = 95.; common[i] = False
            _, audit = qualify(frame, common, design)
            self.assertIn(part + ':early_censoring', audit['gate']['blockers'])
        frame, common, design = synthetic()
        frame.loc[frame.dead.eq(0), 't_exam'] = 96.
        self.assertTrue(qualify(frame, common, design)[1]['gate']['temporal_evaluation_allowed'])

    def test_deaths_at_zero_and96_count_but_deaths_after96_are_binary_nonevents(self):
        frame, common, design = synthetic()
        event = frame.index[frame.cycle.eq(2003) & frame.dead.eq(1)][:3]
        frame.loc[event, 't_exam'] = [0., 96., 97.]
        _, audit = qualify(frame, common, design)
        self.assertEqual(audit['periods']['late']['common']['events8'], 29)
        self.assertEqual(audit['periods']['late']['common']['nonevents8'], 31)
        self.assertTrue(audit['gate']['temporal_evaluation_allowed'])

    def test_empty_common_period_is_explicitly_blocked_and_mask_layout_is_strict(self):
        frame, common, design = synthetic(); common[frame.cycle.eq(2003)] = False
        _, audit = qualify(frame, common, design)
        self.assertIn('late:empty_primary_domain', audit['gate']['blockers'])
        self.assertIsNone(audit['periods']['late']['common']['effective_n'])
        json.dumps(audit, allow_nan=False)
        for bad in (common[:-1], common.astype(int), np.ones((len(frame), 1), dtype=bool)):
            with self.assertRaisesRegex(ValueError, 'common mask'): qualify(frame, bad, design)


class OutcomeFreeSummaryTests(unittest.TestCase):
    def test_feature_and_log2_protein_moments_are_native_weighted_by_period_without_fitting(self):
        frame, common, design = synthetic()
        with patch('research.mortality.b2m31.model.fit_preprocessing', side_effect=AssertionError('unexpected fit')), \
             patch('research.mortality.b2m31.model.fit_ridge', side_effect=AssertionError('unexpected fit')):
            _, audit = qualify(frame, common, design)
        for part, mask, weight_column in [('train', frame.cycle.ne(2003), 'WTSSCB4Y'),
                                          ('late', frame.cycle.eq(2003), 'WTSSCB2Y')]:
            selected = frame.loc[mask]; w = selected[weight_column].to_numpy()
            summary = audit['periods'][part]['common']
            for key, log in [(name, False) for name in c.CORE] + [(name, True) for name in c.PROTEINS]:
                x = selected[key].to_numpy(); x = np.log2(x) if log else x
                mean = float(np.sum(w * x) / np.sum(w))
                sd = float(np.sqrt(np.sum(w * (x - mean) ** 2) / np.sum(w)))
                group = 'outcome_free_log2_protein_moments' if log else 'outcome_free_feature_moments'
                observed = summary[group][key]
                with self.subTest(part=part, key=key):
                    self.assertEqual(observed['n'], len(selected))
                    self.assertAlmostEqual(observed['weighted_mean'], mean, places=12)
                    self.assertAlmostEqual(observed['weighted_SD'], sd, places=12)

    def test_moments_do_not_use_outcomes_and_late_features_do_not_change_early_summaries(self):
        frame, common, design = synthetic(); _, before = qualify(frame, common, design)
        changed = frame.copy(deep=True); changed['dead'] = 1 - changed.dead
        changed['t_exam'] = np.where(changed.dead.eq(1), 48., 177.)
        _, after_outcomes = qualify(changed, common, design)
        groups = ('outcome_free_feature_moments', 'outcome_free_log2_protein_moments')
        for part in ('train', 'late'):
            for group in groups:
                self.assertEqual(before['periods'][part]['common'][group], after_outcomes['periods'][part]['common'][group])
        changed = frame.copy(deep=True); late = changed.cycle.eq(2003)
        changed.loc[late, list(c.CORE)] += 1000.
        changed.loc[late, list(c.PROTEINS)] *= 8.
        _, after_late = qualify(changed, common, design)
        for group in groups:
            self.assertEqual(before['periods']['train']['common'][group], after_late['periods']['train']['common'][group])
        self.assertNotEqual(before['periods']['late']['common']['outcome_free_feature_moments'],
                            after_late['periods']['late']['common']['outcome_free_feature_moments'])


class SourceWrapperTests(unittest.TestCase):
    def test_protocol_schema_hash_and_fingerprint_mismatch_block_before_any_parser(self):
        mutations = [('schema', 'forward32-retrieval/2'), ('plan_commit', 'other-plan'),
                     ('lock_sha256', '0' * 64), ('base_manifest_sha256', '0' * 64),
                     ('source_fingerprints', {}),
                     ('source_fingerprints', {'1999_SSCARD_A.xpt': 'b' * 64})]
        with tempfile.TemporaryDirectory() as directory:
            _, good = wrapper_fixture(directory)
            for key, value in mutations:
                bad = copy.deepcopy(good); bad[key] = value
                (Path(directory) / 'forward32_manifest.json').write_text(json.dumps(bad))
                with self.subTest(key=key), patch.object(c, 'check_study', return_value=LOCK), \
                     patch.object(c.prior, 'load') as parse_prior, \
                     patch.object(c, 'read_xpt_checked') as parse_xpt:
                    with self.assertRaisesRegex(ValueError, 'not bound'): c.load(directory)
                    parse_prior.assert_not_called(); parse_xpt.assert_not_called()

    def test_changed_manifest_bytes_block_before_prior_source_parse(self):
        with tempfile.TemporaryDirectory() as directory:
            wrapper_fixture(directory)
            path = Path(directory) / 'manifest.json'; path.write_bytes(path.read_bytes() + b'\n')
            with patch.object(c, 'check_study', return_value=LOCK), patch.object(c.prior, 'load') as load:
                with self.assertRaisesRegex(ValueError, 'not bound'): c.load(directory)
                load.assert_not_called()

    def test_prior_manifest_disagreement_blocks_before_assay_parse(self):
        with tempfile.TemporaryDirectory() as directory:
            base, _ = wrapper_fixture(directory); frame, common, design = synthetic()
            bad_base = copy.deepcopy(base); bad_base['schema'] = 'different'
            with patch.object(c, 'check_study', return_value=LOCK), \
                 patch.object(c.prior, 'load', return_value=(frame, common, common, {}, bad_base)), \
                 patch.object(c, 'read_xpt_checked') as parse_xpt:
                with self.assertRaisesRegex(ValueError, 'Base source manifest mismatch'): c.load(directory)
                parse_xpt.assert_not_called()

    def test_current_assay_read_hash_is_bound_before_exporting_observed_marginals(self):
        with tempfile.TemporaryDirectory() as directory:
            base, _ = wrapper_fixture(directory); frame, common, design = synthetic()
            fake_assay = pd.DataFrame(); fake_assay.attrs['xpt_zero_audit'] = {'sha256': 'b' * 64}
            with patch.object(c, 'check_study', return_value=LOCK), \
                 patch.object(c.prior, 'load', return_value=(frame, common, common, {'full_DEMO_masked_design': design}, base)), \
                 patch.object(c, 'read_xpt_checked', return_value=fake_assay):
                with self.assertRaisesRegex(ValueError, 'Assay bytes changed'): c.load(directory)

    def test_observed_full_release_flag_and_weight_aggregates_are_exported_from_synthetic_source(self):
        with tempfile.TemporaryDirectory() as directory:
            base, wrapper = wrapper_fixture(directory); frame, common, design = synthetic()
            fake_assay = pd.DataFrame(dict(SSCYSTL=[0., 0., 1., 2., np.nan],
                SSB2ML=[0., 1., 1., 2., np.nan], SSPRIS=[0., 0., 0., 1., np.nan],
                WTSSCB4Y=[10., 0., np.nan, 30., np.nan], WTSSCB2Y=[1., 2., 0., 4., np.nan]))
            fake_assay.attrs['xpt_zero_audit'] = {'sha256': 'a' * 64}
            with patch.object(c, 'check_study', return_value=LOCK), \
                 patch.object(c.prior, 'load', return_value=(frame, common, common, {'full_DEMO_masked_design': design}, base)), \
                 patch.object(c, 'read_xpt_checked', return_value=fake_assay):
                d, actual_common, audit, manifest = c.load(directory)
            observed = audit['observed_full_SSCARD_marginals']
            self.assertEqual(observed['rows'], 5)
            self.assertEqual(observed['flags']['SSCYSTL'], {'0': 2, '1': 1, '2': 1, 'missing': 1})
            self.assertEqual(observed['flags']['SSB2ML'], {'0': 1, '1': 2, '2': 1, 'missing': 1})
            self.assertEqual(observed['pristine'], {'0': 3, '1': 1, 'missing': 1})
            self.assertEqual(observed['weights']['WTSSCB4Y'], {'positive': 2, 'missing': 2, 'zero': 1})
            self.assertEqual(observed['weights']['WTSSCB2Y'], {'positive': 3, 'missing': 1, 'zero': 1})
            self.assertEqual(observed['source_sha256'], 'a' * 64)
            self.assertFalse(observed['individual_batch_or_assay_date_fields_provided'])
            self.assertEqual(manifest, dict(schema='forward32-source-manifest/1', intake=wrapper, source=base))
            self.assertEqual(audit['original31_common_n'], len(frame))
            self.assertTrue(audit['prior31_reported_sensitivity_not_reanalysed'])
            np.testing.assert_array_equal(actual_common, common)
            self.assertEqual(len(d), len(frame))

    def test_retrieval_binds_checked_prior_manifest_without_new_source_requests(self):
        with tempfile.TemporaryDirectory() as directory:
            out = Path(directory) / 'synthetic-intake'; calls = []
            base = dict(schema='synthetic-base-manifest/1', files=[dict(file='synthetic.xpt', sha256='c' * 64)])
            def checked():
                calls.append('study'); return LOCK
            def prior_retrieve(target):
                self.assertEqual(calls, ['study']); calls.append('prior')
                target.mkdir(); (target / 'manifest.json').write_text(json.dumps(base)); return base
            with patch.object(r, 'check_study', side_effect=checked), \
                 patch.object(r.prior, 'retrieve', side_effect=prior_retrieve) as prior:
                wrapper = r.retrieve(out)
            self.assertEqual(calls, ['study', 'prior']); prior.assert_called_once_with(out)
            self.assertEqual(wrapper['plan_commit'], '377778c')
            self.assertEqual(wrapper['source_fingerprints'], {'synthetic.xpt': 'c' * 64})
            self.assertEqual(wrapper['base_manifest_sha256'], r.sha(out / 'manifest.json'))
            self.assertEqual(json.loads((out / 'forward32_manifest.json').read_text()), wrapper)


if __name__ == '__main__':
    unittest.main()
