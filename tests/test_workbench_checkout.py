"""Byte-preserving checkout regression; no clinical observations or fitting."""
from hashlib import sha256
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest
from research.mortality.workbench.runtime import BUNDLE_PATH, BUNDLE_SHA256, load_bundle

ROOT = Path(__file__).resolve().parents[1]


class CheckoutTests(unittest.TestCase):
    def test_checkout_preserves_model_with_autocrlf_enabled(self):
        if not shutil.which('git'):
            self.skipTest('Git is required only for the checkout regression')
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory); repo = root / 'repo'; repo.mkdir()
            target = root / 'checkout'; target.mkdir()
            relative = Path('research/mortality/reanalysis06/model_bundle.json')
            model = repo / relative; model.parent.mkdir(parents=True)
            model.write_bytes(BUNDLE_PATH.read_bytes())
            (repo / '.gitattributes').write_bytes((ROOT / '.gitattributes').read_bytes())
            def git(*args):
                return subprocess.run(['git', *args], cwd=repo, capture_output=True,
                                      text=True, check=True, timeout=20)
            git('init', '--quiet')
            git('config', 'core.autocrlf', 'true')
            git('add', '.gitattributes', relative.as_posix())
            git('checkout-index', '--all', '--force', '--prefix=' + target.as_posix() + '/')
            self.assertEqual(sha256((target / relative).read_bytes()).hexdigest(), BUNDLE_SHA256)
            self.assertIn('text: unset', git('check-attr', 'text', '--', relative.as_posix()).stdout)

    def test_source_policy_preserves_hashed_metadata(self):
        policy = (ROOT / '.gitattributes').read_text(encoding='utf-8')
        for line in ('* text=auto eol=lf', '*.json -text', '*.csv -text'):
            self.assertIn(line, policy)

    def test_line_ending_corruption_still_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            changed = Path(directory) / 'model.json'
            changed.write_bytes(BUNDLE_PATH.read_bytes().replace(b'\n', b'\r\n'))
            with self.assertRaises(ValueError):
                load_bundle(changed)

    def test_CI_native_command_failures_are_fatal(self):
        workflow = (ROOT / '.github/workflows/mortality-workbench.yml').read_text(encoding='utf-8')
        self.assertIn('defaults:\n  run:\n    shell: bash', workflow)
        self.assertIn("paths: ['.gitattributes'", workflow)


if __name__ == '__main__':
    unittest.main()
