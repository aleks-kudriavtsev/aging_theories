"""Execute aggregate-only cause qualification; this stage never fits a model."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import platform
from .retrieve import ROOT, PLAN_COMMIT, check_study, sha, verify_wrapper


def dump(value):
    return (json.dumps(value, ensure_ascii=False, indent=2,
                       allow_nan=False) + '\n').encode('utf-8')


def run(source, out):
    # Lazy imports keep the wrapper/protocol contract independently testable.
    from ..b2m31.cohort import load
    from .gate import qualify

    source, out = Path(source), Path(out)
    if out.exists():
        raise ValueError('New qualification output directory required')
    lock = check_study()
    base, wrapper = verify_wrapper(source)
    eligible, common, _, audit, loaded_base = load(source)
    if loaded_base != base:
        raise ValueError('Source manifest changed between verification and parsing')
    rules = json.loads((ROOT / 'SUPPORT_RULES.json').read_text())
    definitions = json.loads((ROOT / 'CAUSE_DEFINITIONS.json').read_text())
    result = qualify(eligible, common, audit, rules, definitions)
    result.update(schema='cause33-qualified-counts/1', plan_commit=PLAN_COMMIT,
        lock_sha256=sha(ROOT / 'LOCK.json'),
        runtime_python=platform.python_version(),
        prior_cohort_use='all_cycles_previously_used_no_independent_validation',
        comparison_models_fitted=False, predictive_metrics_computed=False,
        independent_binary_probabilities_summed=False,
        joint_integration_and_Fine_Gray='deferred_not_implemented_or_run',
        historical_scientific_files_verified=len(lock['historical_scientific_python_sha256']),
        scientific_python_sha256={str(p.relative_to(ROOT.parent)): sha(p)
            for p in sorted(ROOT.parent.rglob('*.py'))})
    manifest = dict(schema='cause33-qualified-sources/1', wrapper=wrapper,
        base_manifest=base, inherited_source_domain_audit=audit,
        retained_endpoint_contract_sha256=sha(ROOT.parent / 'endpoint_contract11.py'),
        endpoint_review_sha256=sha(ROOT / 'ENDPOINT_REVIEW.md'),
        dictionary_status='official_indexed_text_and_retained_review_not_new_PDF_bytes',
        participant_rows_exported=False, ZIP_downloaded=False)
    outputs = {'GATE_RESULTS.json': result, 'SOURCE_MANIFEST.json': manifest}
    # Serialize completely before writing: NaN or unsafe objects cannot leave
    # a partial success file that could be mistaken for completed qualification.
    encoded = {name: dump(value) for name, value in outputs.items()}
    out.mkdir(parents=True)
    for name, content in encoded.items():
        (out / name).write_bytes(content)
        print('CAUSE33_STDOUT_FILE_BEGIN ' + name, flush=True)
        print(content.decode(), end='', flush=True)
        print('CAUSE33_STDOUT_FILE_END ' + name, flush=True)
    print('CAUSE33_GATE_COMPLETE no_comparative_fits', flush=True)
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    arguments = parser.parse_args()
    run(arguments.source, arguments.out)
