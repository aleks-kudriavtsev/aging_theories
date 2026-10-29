"""No production simulations: runner contracts use explicit mocks only."""
import contextlib
import io
import json
from pathlib import Path
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

    def test_all_frozen_case_inputs_passed_unchanged_and_exports_are_exact(self):
        with tempfile.TemporaryDirectory() as directory:
            out = Path(directory) / 'new_results'
            qualification = {'source_ready_candidate_ids': []}
            stdout = io.StringIO()
            with patch.object(a.precision, 'plan', return_value={'synthetic': True}) as planner, \
                    patch.object(a.precision, 'simulate_case', side_effect=toy) as simulate, \
                    patch.object(a.gate, 'qualify_all', return_value=qualification) as qualify, \
                    contextlib.redirect_stdout(stdout):
                result = a.run(out)
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
