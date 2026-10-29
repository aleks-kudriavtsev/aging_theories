"""Bind the preserved official public intake to the cause33A count protocol."""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
from ..b2m31 import retrieve as prior

ROOT = Path(__file__).parent
PLAN_COMMIT = '82eacda12eced3c0105019fee09236b647821f45'
LOCAL_METADATA = frozenset({'GATE_PLAN.md', 'CAUSE_DEFINITIONS.json',
    'SUPPORT_RULES.json', 'ENDPOINT_REVIEW.md', 'DESIGN_REVIEW.md'})


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def check_study():
    lock = json.loads((ROOT / 'LOCK.json').read_text())
    if lock.get('schema') != 'cause33-gate-lock/1':
        raise ValueError('Wrong cause33 protocol lock')
    for name, digest in lock['metadata_sha256'].items():
        path = ROOT / name if name in LOCAL_METADATA else ROOT.parent / name
        if sha(path) != digest:
            raise ValueError('Cause33 locked metadata changed')
    for name, digest in lock['historical_scientific_python_sha256'].items():
        if sha(ROOT.parent / name) != digest:
            raise ValueError('Previous scientific code changed')
    return lock


def retrieve(out):
    check_study()
    out = Path(out)
    base = prior.retrieve(out)
    wrapper = dict(schema='cause33-retrieval/1', plan_commit=PLAN_COMMIT,
        lock_sha256=sha(ROOT / 'LOCK.json'),
        base_manifest_sha256=sha(out / 'manifest.json'),
        source_fingerprints={row['file']: row['sha256'] for row in base['files']})
    (out / 'cause33_manifest.json').write_text(
        json.dumps(wrapper, indent=2, allow_nan=False) + '\n')
    return wrapper


def verify_wrapper(source):
    check_study()
    source = Path(source)
    base = json.loads((source / 'manifest.json').read_text())
    wrapper = json.loads((source / 'cause33_manifest.json').read_text())
    if (wrapper.get('schema') != 'cause33-retrieval/1'
            or wrapper.get('plan_commit') != PLAN_COMMIT
            or wrapper.get('lock_sha256') != sha(ROOT / 'LOCK.json')
            or wrapper.get('base_manifest_sha256') != sha(source / 'manifest.json')
            or wrapper.get('source_fingerprints') != {
                row['file']: row['sha256'] for row in base['files']}):
        raise ValueError('Source wrapper not bound to cause33 protocol')
    return base, wrapper


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    retrieve(parser.parse_args().out)
