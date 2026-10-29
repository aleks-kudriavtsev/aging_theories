"""Independent publication-auditor tests; no producer imports or simulations."""
import copy
import importlib.util
import json
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

REPO = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location(
    'replication34_publication_auditor', REPO / 'tools/verify_replication34_publication.py')
a = importlib.util.module_from_spec(spec)
spec.loader.exec_module(a)
LOCAL = REPO / 'research/mortality/replication34'


class PublicationAuditTests(unittest.TestCase):
    def summary(self):
        return a.read_json(LOCAL / 'RESULTS_SUMMARY.json')

    def test_only_explicit_version_and_runtime_hash_changes_allowed(self):
        original = self.summary(); current = copy.deepcopy(original)
        current['versions']['python'] = '3.13.7'
        current['scientific_python_sha256']['replication34/analyze.py'] = 'a' * 64
        self.assertEqual(a.summary_differences(original, current), [])
        current['qualified_algorithm_source_ids'] = ['invented_ready_cohort']
        self.assertEqual(a.summary_differences(original, current), ['/qualified_algorithm_source_ids/length'])

    def test_missing_keys_extra_keys_and_numeric_type_substitution_rejected(self):
        original = self.summary()
        for change in ('remove', 'add', 'type'):
            current = copy.deepcopy(original)
            if change == 'remove':
                del current['model_fits_performed']
            elif change == 'add':
                current['new_claim'] = True
            else:
                current['model_fits_performed'] = 0
            self.assertTrue(a.summary_differences(original, current), change)

    def test_duplicate_or_nonfinite_json_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'bad.json'
            for text in ('{"x": 1, "x": 2}', '{"x": NaN}', '{"x": Infinity}'):
                path.write_text(text)
                with self.assertRaises(ValueError):
                    a.read_json(path)

    def test_locked_path_escape_and_symlink_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory); (root / 'safe').write_text('safe')
            (root / 'linked').symlink_to(root / 'safe')
            for name in ('../safe', '/tmp/safe', 'missing', 'linked'):
                with self.assertRaises(ValueError, msg=name):
                    a.checked_path(root, name)
            self.assertEqual(a.checked_path(root, 'safe'), root / 'safe')

    def test_source_sibling_references_route_to_history_not_main(self):
        with tempfile.TemporaryDirectory() as directory:
            base = Path(directory); local = base / 'publication'; history = base / 'history'
            local.mkdir(); (history / 'freeze30').mkdir(parents=True)
            (local / 'review.md').write_text('current')
            (history / 'freeze30' / 'review.md').write_text('pinned historical')
            self.assertEqual(a.source_path(local, history, 'review.md'), local / 'review.md')
            self.assertEqual(a.source_path(local, history, '../freeze30/review.md'), history / 'freeze30/review.md')
            with self.assertRaises(ValueError):
                a.source_path(local, history, '../../escape')

    def test_original_independent_math_file_is_pinned_before_import(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory); (root / 'tools').mkdir()
            (root / 'tools/verify_replication34_aggregates.py').write_text(
                "raise AssertionError('must not execute an altered auditor')\n")
            with self.assertRaisesRegex(ValueError, 'auditor changed'):
                a.load_independent_math(root)

    def test_marked_stdout_only_strips_transport_prefix_and_requires_unique_block(self):
        independent = a.load_independent_math(REPO)
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'log'
            payload = '{"text": "\ufeffinside payload remains"}\n'
            block = 'REP34_FILE_BEGIN demo.json\n' + payload + 'REP34_FILE_END demo.json\n'
            path.write_text(''.join('\ufeff2026-10-29T18:15:16.0632843Z ' + line
                                    for line in block.splitlines(keepends=True)))
            self.assertEqual(independent.marked(independent.log_text(path), 'demo.json'), payload)
            with self.assertRaisesRegex(ValueError, 'Nonunique'):
                independent.marked(block + block, 'demo.json')

    def test_immutable_archive_hashes_and_independent_math(self):
        independent = a.load_independent_math(REPO)
        self.assertTrue(all(a.sha(LOCAL / name) == digest for name, digest in a.ORIGINAL_SHA.items()))
        data = {name: a.read_json(LOCAL / name) for name in a.ORIGINAL_SHA}
        audit = a.Audit(); a.audit_math(audit, independent, LOCAL, data)
        self.assertTrue(all(row['passed'] for row in audit.checks),
                        [row for row in audit.checks if not row['passed']])
        altered = copy.deepcopy(data)
        altered['PRECISION_RESULTS.json']['cases'][0]['result']['conditional_n_for_half_width'] += 1
        audit = a.Audit(); a.audit_math(audit, independent, LOCAL, altered)
        self.assertFalse(all(row['passed'] for row in audit.checks))

    @unittest.skipUnless(os.environ.get('REP34_HISTORY_ROOT'), 'requires explicit pinned history checkout')
    def test_full_publication_summary_and_perturbed_payload_fail_closed(self):
        history = Path(os.environ['REP34_HISTORY_ROOT'])
        # Copy immutable aggregate fixtures; only code/version metadata is adapted.
        # This is an audit test fixture, never a new production execution claim.
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)
            for name in a.ORIGINAL_SHA:
                (output / name).write_bytes((LOCAL / name).read_bytes())
            summary = self.summary()
            for name in a.CURRENT_PATHS:
                summary['scientific_python_sha256'][name] = a.sha(REPO / 'research/mortality' / name)
            (output / 'RESULTS_SUMMARY.json').write_text(json.dumps(summary))
            report = a.verify(REPO, history, output)
            self.assertTrue(report['all_passed'], report['failed'])
            payload = output / 'PRECISION_RESULTS.json'
            payload.write_bytes(payload.read_bytes() + b' ')
            report = a.verify(REPO, history, output)
            self.assertFalse(report['all_passed'])
            self.assertIn('byte-identical reproduced scientific payload PRECISION_RESULTS.json', report['failed'])
            with patch.object(a, 'git_head', return_value='0' * 40):
                with self.assertRaisesRegex(ValueError, 'Historical checkout must be exactly'):
                    a.verify(REPO, history, output)


if __name__ == '__main__':
    unittest.main()
