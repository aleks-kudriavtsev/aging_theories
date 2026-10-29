"""Synthetic protocol binding; no source downloads or participant data."""
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
from research.mortality.cause33 import retrieve as r


class IntakeTests(unittest.TestCase):
    def fixture(self, directory):
        directory = Path(directory)
        base = dict(files=[dict(file='synthetic.xpt', sha256='a' * 64)])
        (directory / 'manifest.json').write_text(json.dumps(base))
        wrapper = dict(schema='cause33-retrieval/1', plan_commit=r.PLAN_COMMIT,
            lock_sha256=r.sha(r.ROOT / 'LOCK.json'),
            base_manifest_sha256=r.sha(directory / 'manifest.json'),
            source_fingerprints={'synthetic.xpt': 'a' * 64})
        (directory / 'cause33_manifest.json').write_text(json.dumps(wrapper))
        return base, wrapper

    def test_actual_plan_metadata_and_all95_old_modules_remain_locked(self):
        lock = r.check_study()
        self.assertEqual(len(lock['historical_scientific_python_sha256']), 95)
        self.assertEqual(len(lock['metadata_sha256']), 18)

    def test_wrapper_verifies_exact_manifest_plan_and_fingerprints(self):
        with tempfile.TemporaryDirectory() as directory:
            expected = self.fixture(directory)
            self.assertEqual(r.verify_wrapper(directory), expected)

    def test_every_wrapper_contract_mutation_is_rejected(self):
        for field, value in [('schema', 'wrong'), ('plan_commit', 'wrong'),
                ('lock_sha256', 'b' * 64), ('base_manifest_sha256', 'b' * 64),
                ('source_fingerprints', {})]:
            with self.subTest(field=field), tempfile.TemporaryDirectory() as directory:
                _, wrapper = self.fixture(directory)
                wrapper[field] = value
                (Path(directory) / 'cause33_manifest.json').write_text(json.dumps(wrapper))
                with self.assertRaisesRegex(ValueError, 'not bound'):
                    r.verify_wrapper(directory)

    def test_changed_base_manifest_is_not_accepted(self):
        with tempfile.TemporaryDirectory() as directory:
            self.fixture(directory)
            (Path(directory) / 'manifest.json').write_text('{"files": []}')
            with self.assertRaisesRegex(ValueError, 'not bound'):
                r.verify_wrapper(directory)

    def test_wrapper_retrieval_only_calls_preserved_intake(self):
        with tempfile.TemporaryDirectory() as directory:
            base, _ = self.fixture(directory)
            with patch.object(r.prior, 'retrieve', return_value=base) as retrieve:
                wrapper = r.retrieve(directory)
            retrieve.assert_called_once_with(Path(directory))
            self.assertEqual(wrapper['plan_commit'], r.PLAN_COMMIT)
            self.assertEqual(r.verify_wrapper(directory)[1], wrapper)


if __name__ == '__main__':
    unittest.main()
