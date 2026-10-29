"""Runner isolation and complete aggregate export, using synthetic mocks only."""
import contextlib
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import pandas as pd
from research.mortality.cause33 import qualify as runner
from research.mortality.b2m31 import cohort
from research.mortality.cause33 import gate


class RunnerTests(unittest.TestCase):
    def test_export_contains_aggregates_and_source_audit_without_input_frame(self):
        frame = pd.DataFrame({'SEQN': [987654321]})
        base = {'files': []}
        result = {'gate': {'later_comparisons_allowed': False}, 'proceed_codes': []}
        with tempfile.TemporaryDirectory() as directory:
            source, out = Path(directory), Path(directory) / 'new_outputs'
            with patch.object(runner, 'verify_wrapper', return_value=(base, {})), \
                    patch.object(cohort, 'load', return_value=(frame, [True], None, {}, base)), \
                    patch.object(gate, 'qualify', return_value=result) as qualify, \
                    contextlib.redirect_stdout(io.StringIO()):
                runner.run(source, out)
            self.assertIs(qualify.call_args.args[0], frame)
            self.assertEqual({p.name for p in out.iterdir()}, {'GATE_RESULTS.json', 'SOURCE_MANIFEST.json'})
            exported = '\n'.join(p.read_text() for p in out.iterdir())
            self.assertNotIn('987654321', exported)
            summary = json.loads((out / 'GATE_RESULTS.json').read_text())
            self.assertFalse(summary['comparison_models_fitted'])
            self.assertFalse(summary['independent_binary_probabilities_summed'])
            self.assertFalse(summary['gate']['later_comparisons_allowed'])

    def test_existing_output_directory_rejected_before_source_load(self):
        with tempfile.TemporaryDirectory() as directory, \
                patch.object(cohort, 'load', side_effect=AssertionError('must not load')):
            with self.assertRaisesRegex(ValueError, 'New qualification'):
                runner.run(directory, directory)

    def test_manifest_change_during_parsing_leaves_no_output(self):
        with tempfile.TemporaryDirectory() as directory:
            out = Path(directory) / 'new_outputs'
            with patch.object(runner, 'verify_wrapper', return_value=({'files': []}, {})), \
                    patch.object(cohort, 'load', return_value=(None, None, None, {}, {'files': ['changed']})):
                with self.assertRaisesRegex(ValueError, 'manifest changed'):
                    runner.run(directory, out)
            self.assertFalse(out.exists())

    def test_nonfinite_result_does_not_leave_partial_success_exports(self):
        with tempfile.TemporaryDirectory() as directory:
            out = Path(directory) / 'new_outputs'
            base = {'files': []}
            with patch.object(runner, 'verify_wrapper', return_value=(base, {})), \
                    patch.object(cohort, 'load', return_value=(None, None, None, {}, base)), \
                    patch.object(gate, 'qualify', return_value={'invalid': float('nan')}):
                with self.assertRaises(ValueError):
                    runner.run(directory, out)
            self.assertFalse(out.exists())


if __name__ == '__main__':
    unittest.main()
