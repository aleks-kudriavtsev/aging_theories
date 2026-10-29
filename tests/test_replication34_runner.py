"""No production simulations: runner contracts use explicit mocks only."""
import contextlib
import io
import json
import os
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch
from research.mortality.replication34 import analyze as a


def toy(**kwargs):
    return dict(MC_mean_d=0., exact_mean_d=0., MCSE_mean_d=0.,
        MC_variance_mean_d=0., exact_variance_mean_d=0.,
        approximate_fourth_moment_MCSE_variance=0.,
        MC_mean_event_count=30., exact_expected_event_count=30.,
        MCSE_mean_event_count=1., replicates=kwargs['replicates'])


class RunnerTests(unittest.TestCase):
    def test_protocol_source_locks_and103_old_scientific_modules(self):
        lock = a.check_study()
        self.assertEqual(len(lock['metadata_sha256']), 22)
        self.assertEqual(len(lock['historical_scientific_python_sha256']), 103)
        self.assertEqual(lock['base_commit'], a.HISTORY_COMMIT)

    def test_explicit_and_environment_history_roots_and_legacy_default(self):
        history = a._history_root()
        with patch.dict(os.environ, {'REP34_HISTORY_ROOT': str(history)}):
            self.assertEqual(a._history_root(), history)
        with patch.dict(os.environ, {'REP34_HISTORY_ROOT': '/missing/ignored_history'}):
            self.assertEqual(a._history_root(history), history)
        # The old single-checkout layout still defaults to ROOT.parents[2].
        with patch.dict(os.environ, {}, clear=True), \
                patch.object(a, 'ROOT', history / 'research' / 'mortality' / 'replication34'):
            self.assertEqual(a._history_root(), history)

    def test_missing_empty_nonrepository_and_subdirectory_history_roots_fail(self):
        history = a._history_root()
        with tempfile.TemporaryDirectory() as directory:
            with self.assertRaisesRegex(ValueError, 'does not exist'):
                a.check_study(Path(directory) / 'missing')
            with self.assertRaisesRegex(ValueError, 'cannot be verified'):
                a.check_study(directory)
        with self.assertRaisesRegex(ValueError, 'repository root required'):
            a.check_study('')
        with self.assertRaisesRegex(ValueError, 'must be the repository root'):
            a.check_study(history / 'research')

    def test_wrong_history_HEAD_fails_before_any_planning_or_simulation(self):
        history = a._history_root()
        responses = [subprocess.CompletedProcess([], 0, str(history) + '\n'),
                     subprocess.CompletedProcess([], 0, '0' * 40 + '\n')]
        with tempfile.TemporaryDirectory() as directory, \
                patch.object(a.subprocess, 'run', side_effect=responses), \
                patch.object(a.precision, 'plan') as planner, \
                patch.object(a.precision, 'simulate_case') as simulate:
            out = Path(directory) / 'new_results'
            with self.assertRaisesRegex(ValueError, 'HEAD differs'):
                a.run(out, history_root=history)
            self.assertFalse(out.exists())
            planner.assert_not_called()
            simulate.assert_not_called()

    def test_all_pins_resolve_to_their_designated_current_or_historical_root(self):
        history = a._history_root()
        historical_mortality = history / 'research' / 'mortality'
        with patch.object(a, 'sha', wraps=a.sha) as checked:
            lock = a.check_study(history)
        paths = [Path(call.args[0]).resolve() for call in checked.call_args_list]
        local = set(lock['local_metadata'])
        for name in lock['metadata_sha256']:
            expected = a.ROOT / name if name in local else historical_mortality / name
            self.assertIn(expected.resolve(), paths)
        for name in lock['historical_scientific_python_sha256']:
            self.assertIn((historical_mortality / name).resolve(), paths)
        manifest = json.loads((a.ROOT / 'SOURCE_MANIFEST.json').read_text())
        self.assertEqual(len(manifest['review_sha256']), 7)
        for name in manifest['review_sha256']:
            expected = (historical_mortality / name[3:] if name.startswith('../')
                        else a.ROOT / name)
            self.assertIn(expected.resolve(), paths)
        catalog = historical_mortality / 'freeze30' / 'CANDIDATE_MANIFEST.json'
        # Both the protocol pin and source-catalog pin are checked.
        self.assertEqual(paths.count(catalog.resolve()), 2)
        self.assertEqual(len(paths), 22 + 103 + 2 + 7 + 1)

    def test_changed_local_metadata_historical_metadata_code_review_and_catalog_fail(self):
        history = a._history_root()
        historical_mortality = history / 'research' / 'mortality'
        original_sha = a.sha
        cases = [
            (a.ROOT / 'SCENARIOS.json', 1, 'metadata changed'),
            (historical_mortality / 'cause33b' / 'RESULTS_SUMMARY.json', 1, 'metadata changed'),
            (historical_mortality / 'engine.py', 1, 'scientific code changed'),
            (a.ROOT / 'PRE_RESULT_SOURCE_ADDENDUM.md', 1, 'review evidence changed'),
            (historical_mortality / 'freeze30' / 'QUALIFICATION_PUBLIC_RU.md', 2,
             'review evidence changed'),
            (historical_mortality / 'freeze30' / 'CANDIDATE_MANIFEST.json', 2,
             'source catalogue changed'),
        ]
        for target, failure_occurrence, message in cases:
            with self.subTest(path=target, failure_occurrence=failure_occurrence):
                occurrences = 0

                def altered_sha(path):
                    nonlocal occurrences
                    if Path(path).resolve() == target.resolve():
                        occurrences += 1
                        if occurrences == failure_occurrence:
                            return '0' * 64
                    return original_sha(path)

                with patch.object(a, 'sha', side_effect=altered_sha):
                    with self.assertRaisesRegex(ValueError, message):
                        a.check_study(history)

    def test_review_path_cannot_escape_current_or_historical_root(self):
        history = a._history_root()
        for name in ('../../outside.md', '/outside.md'):
            with self.subTest(name=name), self.assertRaisesRegex(ValueError, 'escapes'):
                a._source_review_path(name, history)

    def test_all_frozen_case_inputs_passed_unchanged_and_exports_are_exact(self):
        with tempfile.TemporaryDirectory() as directory:
            out = Path(directory) / 'new_results'
            qualification = {'source_ready_candidate_ids': []}
            stdout = io.StringIO()
            with patch.object(a.precision, 'plan', return_value={'synthetic': True}) as planner, \
                    patch.object(a.precision, 'simulate_case', side_effect=toy) as simulate, \
                    patch.object(a.gate, 'qualify_all', return_value=qualification) as qualify, \
                    patch.object(Path, 'rglob', side_effect=AssertionError('No main-tree scan')), \
                    contextlib.redirect_stdout(stdout):
                result = a.run(out, history_root=a._history_root())
            self.assertEqual(planner.call_count, 12)
            self.assertEqual(simulate.call_count, 4)
            self.assertEqual([call.kwargs['seed'] for call in simulate.call_args_list],
                             [20261120, 20261121, 20261122, 20261123])
            self.assertTrue(all(call.kwargs['replicates'] == 20000 for call in simulate.call_args_list))
            self.assertEqual([call.kwargs['mode'] for call in qualify.call_args_list],
                             ['algorithm_replication', 'locked_source_transport'])
            self.assertEqual(len(list(out.iterdir())), 5)
            self.assertFalse(result['model_fits_performed'])
            self.assertFalse(result['statistical_adequacy_for_refitted_algorithm_established'])
            lock = json.loads((a.ROOT / 'LOCK.json').read_text())
            expected_science = dict(lock['historical_scientific_python_sha256'])
            expected_science.update({'replication34/' + name: a.sha(a.ROOT / name)
                                     for name in a.CURRENT_SCIENTIFIC_FILES})
            self.assertEqual(result['scientific_python_sha256'], expected_science)
            self.assertEqual(len(result['scientific_python_sha256']), 107)
            for path in out.iterdir():
                text = path.read_text()
                self.assertIn('REP34_FILE_BEGIN ' + path.name + '\n' + text + 'REP34_FILE_END ' + path.name,
                              stdout.getvalue())
                if path.name != 'RESULTS_SUMMARY.json':
                    self.assertEqual(a.sha(path), result['outputs_sha256'][path.name])

    def test_existing_output_stops_before_metadata_or_simulation(self):
        with tempfile.TemporaryDirectory() as directory, \
                patch.object(a, 'check_study', side_effect=AssertionError('must not run')):
            with self.assertRaisesRegex(ValueError, 'New aggregate output'):
                a.run(directory)

    def test_nonfinite_mock_result_leaves_no_partial_output(self):
        with tempfile.TemporaryDirectory() as directory:
            out = Path(directory) / 'new_results'
            with patch.object(a.precision, 'plan', return_value={'bad': float('nan')}), \
                    patch.object(a.precision, 'simulate_case', side_effect=toy), \
                    patch.object(a.gate, 'qualify_all', return_value={'source_ready_candidate_ids': []}):
                with self.assertRaises(ValueError):
                    a.run(out)
            self.assertFalse(out.exists())

    def test_MC_agreement_does_not_assert_nominal_coverage_or_hide_disagreement(self):
        self.assertTrue(a._MC_check(1., 1., 0.)['predeclared_six_MCSE_numerical_check_passed'])
        self.assertFalse(a._MC_check(1.1, 1., 0.)['predeclared_six_MCSE_numerical_check_passed'])
        self.assertFalse(a._MC_check(.1, 0., .001)['predeclared_six_MCSE_numerical_check_passed'])


if __name__ == '__main__':
    unittest.main()
