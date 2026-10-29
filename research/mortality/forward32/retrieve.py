"""Bind the preserved public31 intake to the new forward-calendar protocol."""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
from ..b2m31 import retrieve as prior

ROOT=Path(__file__).parent
PLAN_COMMIT='377778c'


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def check_study():
    lock=json.loads((ROOT/'LOCK.json').read_text())
    for name,digest in lock['metadata_sha256'].items():
        path=ROOT.parent/name if '/' in name else ROOT/name
        if sha(path)!=digest:raise ValueError('Forward32 locked metadata changed')
    for name,digest in lock['historical_scientific_python_sha256'].items():
        if sha(ROOT.parent/name)!=digest:raise ValueError('Previous scientific code changed')
    return lock


def retrieve(out):
    check_study();out=Path(out);base=prior.retrieve(out)
    wrapper=dict(schema='forward32-retrieval/1',plan_commit=PLAN_COMMIT,
        lock_sha256=sha(ROOT/'LOCK.json'),base_manifest_sha256=sha(out/'manifest.json'),
        source_fingerprints={r['file']:r['sha256'] for r in base['files']})
    (out/'forward32_manifest.json').write_text(json.dumps(wrapper,indent=2,allow_nan=False)+'\n')
    return wrapper


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--out',type=Path,required=True)
    retrieve(p.parse_args().out)
