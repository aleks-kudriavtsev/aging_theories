"""Source-wrapper binding without network or participant data."""
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
from research.mortality.cause33b import retrieve as r


class IntakeTests(unittest.TestCase):
    def fixture(self, directory):
        directory = Path(directory)
        base = dict(files=[dict(file='synthetic.xpt', sha256='a' * 64)])
        (directory / 'manifest.json').write_text(json.dumps(base))
        wrapper = dict(schema='cause33b-retrieval/1', plan_commit=r.PLAN_COMMIT,
            lock_sha256=r.sha(r.ROOT / 'LOCK.json'),
            base_manifest_sha256=r.sha(directory / 'manifest.json'),
            source_fingerprints={'synthetic.xpt': 'a' * 64})
        (directory / 'cause33b_manifest.json').write_text(json.dumps(wrapper))
        return base, wrapper

    def test_all20_metadata_and99_old_scientific_files_remain_locked(self):
        lock = r.check_study()
        self.assertEqual(len(lock['metadata_sha256']), 20)
        self.assertEqual(len(lock['historical_scientific_python_sha256']), 99)

    def test_exact_binding_and_every_mutated_contract_field(self):
        with tempfile.TemporaryDirectory() as directory:
            expected = self.fixture(directory)
            self.assertEqual(expected, r.verify_wrapper(directory))
        for key, value in [('schema', 'wrong'), ('plan_commit', 'wrong'),
                ('lock_sha256', 'b' * 64), ('base_manifest_sha256', 'b' * 64),
                ('source_fingerprints', {})]:
            with self.subTest(key=key), tempfile.TemporaryDirectory() as directory:
                _, wrapper = self.fixture(directory); wrapper[key] = value
                (Path(directory) / 'cause33b_manifest.json').write_text(json.dumps(wrapper))
                with self.assertRaisesRegex(ValueError, 'not bound'):
                    r.verify_wrapper(directory)

    def test_retrieval_only_uses_preserved_intake_and_records_exact_base_bytes(self):
        with tempfile.TemporaryDirectory() as directory:
            base, _ = self.fixture(directory)
            with patch.object(r.prior, 'retrieve', return_value=base) as retrieve:
                wrapper = r.retrieve(directory)
            retrieve.assert_called_once_with(Path(directory))
            self.assertEqual(wrapper, r.verify_wrapper(directory)[1])
            (Path(directory) / 'manifest.json').write_text('{"files": []}')
            with self.assertRaisesRegex(ValueError, 'not bound'):
                r.verify_wrapper(directory)


if __name__ == '__main__':
    unittest.main()
