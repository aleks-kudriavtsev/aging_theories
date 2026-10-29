"""Export metadata qualification and stipulated precision/simulation aggregates."""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import platform
import numpy as np
import scipy
from . import gate, precision

ROOT = Path(__file__).parent
PLAN_COMMIT = '4df636ea18ca8bb29b05788c37800438566426ec'


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def encoded(value):
    return (json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + '\n').encode('utf-8')


def check_study():
    lock = json.loads((ROOT / 'LOCK.json').read_text())
    if lock.get('schema') != 'replication34-lock/1':
        raise ValueError('Wrong replication34 protocol lock')
    local = set(lock['local_metadata'])
    for name, digest in lock['metadata_sha256'].items():
        path = ROOT / name if name in local else ROOT.parent / name
        if sha(path) != digest:
            raise ValueError('Locked replication34 metadata changed')
    for name, digest in lock['historical_scientific_python_sha256'].items():
        if sha(ROOT.parent / name) != digest:
            raise ValueError('Prior scientific code changed')
    source_lock = json.loads((ROOT / 'SOURCE_LOCK.json').read_text())
    if (source_lock.get('schema') != 'replication34-source-lock/1'
            or set(source_lock['files_sha256']) != {'CANDIDATES.json', 'SOURCE_MANIFEST.json'}):
        raise ValueError('Wrong source metadata lock')
    for name, digest in source_lock['files_sha256'].items():
        if sha(ROOT / name) != digest:
            raise ValueError('Reviewed source metadata changed')
    manifest = json.loads((ROOT / 'SOURCE_MANIFEST.json').read_text())
    for name, digest in manifest['review_sha256'].items():
        if sha(ROOT / name) != digest:
            raise ValueError('Source review evidence changed')
    return lock


def _MC_check(observed, expected, mcse):
    # This is a numerical simulation check, never a model-performance test.
    if mcse == 0:
        passed = observed == expected
        standardized = None
    else:
        standardized = abs(observed - expected) / mcse
        passed = standardized <= 6.
    return dict(observed=observed, expected=expected, MCSE=mcse,
                absolute_discrepancy_in_MCSE=standardized,
                predeclared_six_MCSE_numerical_check_passed=passed)


def run(out):
    out = Path(out)
    if out.exists():
        raise ValueError('New aggregate output directory required')
    lock = check_study()
    scenario_file = json.loads((ROOT / 'SCENARIOS.json').read_text())
    sim_file = json.loads((ROOT / 'SIMULATION_CASES.json').read_text())
    candidates = json.loads((ROOT / 'CANDIDATES.json').read_text())
    manifest = json.loads((ROOT / 'SOURCE_MANIFEST.json').read_text())
    if (scenario_file.get('schema') != 'replication34-scenarios/1'
            or scenario_file.get('not_estimated_from_33B_effects_or_SEs') is not True
            or sim_file.get('schema') != 'replication34-simulation/1'
            or sim_file.get('predictions_defined_before_outcomes') is not True):
        raise ValueError('Locked assumption/simulation contract required')
    plans = []
    ids = set()
    for row in scenario_file['scenarios']:
        if row['id'] in ids:
            raise ValueError('Duplicate precision scenario')
        ids.add(row['id'])
        plans.append(dict(id=row['id'], assumptions=row, result=precision.plan(row)))
    simulations = []
    ids = set()
    for row in sim_file['cases']:
        if row['id'] in ids:
            raise ValueError('Duplicate simulation case')
        ids.add(row['id'])
        arguments = {k: v for k, v in row.items() if k != 'id'}
        result = precision.simulate_case(**arguments)
        checks = dict(
            mean_D=_MC_check(result['MC_mean_d'], result['exact_mean_d'], result['MCSE_mean_d']),
            variance_Dbar=_MC_check(result['MC_variance_mean_d'], result['exact_variance_mean_d'],
                                   result['approximate_fourth_moment_MCSE_variance']),
            mean_event_count=_MC_check(result['MC_mean_event_count'], result['exact_expected_event_count'],
                                       result['MCSE_mean_event_count']))
        simulations.append(dict(id=row['id'], assumptions=row, result=result,
            moment_checks=checks,
            coverage_95pct_is_reported_not_required_or_tuned=True))
    readiness = dict(schema='replication34-readiness/1',
        algorithm_replication=gate.qualify_all(candidates, manifest, mode='algorithm_replication'),
        locked_source_transport=gate.qualify_all(candidates, manifest, mode='locked_source_transport'),
        source_qualification_is_not_statistical_adequacy=True,
        external_predictions_or_outcome_comparisons_performed=False)
    precision_report = dict(schema='replication34-planning-results/1',
        input_metadata={k: v for k, v in scenario_file.items() if k != 'scenarios'},
        cases=plans, estimates_are_stipulated_scenarios_not_empirical_findings=True,
        adequate_sample_size_or_power_established=False)
    simulation_report = dict(schema='replication34-numerical-validation/1',
        input_metadata={k: v for k, v in sim_file.items() if k != 'cases'},
        cases=simulations,
        all_predeclared_moment_checks_passed=all(
            check['predeclared_six_MCSE_numerical_check_passed']
            for row in simulations for check in row['moment_checks'].values()),
        fitted_algorithm_precision_or_external_validation_performed=False)
    summary = dict(schema='replication34-summary/1', plan_commit=PLAN_COMMIT,
        phase='replication_readiness_and_conditional_precision_only',
        candidates_reviewed=len(candidates),
        qualified_algorithm_source_ids=readiness['algorithm_replication']['source_ready_candidate_ids'],
        qualified_transport_source_ids=readiness['locked_source_transport']['source_ready_candidate_ids'],
        illustrative_precision_scenarios=len(plans), synthetic_validation_cases=len(simulations),
        MonteCarlo_replicates={s['id']: s['result']['replicates'] for s in simulations},
        moment_checks_passed=simulation_report['all_predeclared_moment_checks_passed'],
        before_new_outcomes_named_release_and_precision_goals_required=True,
        statistical_adequacy_for_refitted_algorithm_established=False,
        guaranteed_recruitment_or_clinical_threshold_provided=False,
        new_person_level_data_or_outcomes_read=False, model_fits_performed=False,
        split_search_performed=False, old33B_result_changed=False,
        lab_report_v3_or_presentation_touched=False,
        joint_integration_and_Fine_Gray='deferred_not_implemented_or_run',
        prior_scientific_files_verified=len(lock['historical_scientific_python_sha256']),
        scientific_python_sha256={str(p.relative_to(ROOT.parent)): sha(p)
            for p in sorted(ROOT.parent.rglob('*.py'))},
        metadata_lock_sha256=sha(ROOT / 'LOCK.json'), source_lock_sha256=sha(ROOT / 'SOURCE_LOCK.json'),
        source_manifest_sha256=sha(ROOT / 'SOURCE_MANIFEST.json'),
        versions=dict(python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__))
    payloads = {'PRECISION_RESULTS.json': encoded(precision_report),
        'SIMULATION_VALIDATION.json': encoded(simulation_report),
        'READINESS_RESULTS.json': encoded(readiness),
        'SOURCE_MANIFEST.json': (ROOT / 'SOURCE_MANIFEST.json').read_bytes()}
    summary['outputs_sha256'] = {name: hashlib.sha256(data).hexdigest() for name, data in payloads.items()}
    payloads['RESULTS_SUMMARY.json'] = encoded(summary)
    # No partial success files from nonfinite results or a changed protocol.
    check_study()
    out.mkdir(parents=True)
    for name, data in payloads.items():
        (out / name).write_bytes(data)
        print('REP34_FILE_BEGIN ' + name, flush=True)
        print(data.decode('utf-8'), end='', flush=True)
        print('REP34_FILE_END ' + name, flush=True)
    return summary


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    result = run(parser.parse_args().out)
    if not result['moment_checks_passed']:
        raise SystemExit('Predeclared toy-law numerical moment checks failed; inspect exported diagnostics')
