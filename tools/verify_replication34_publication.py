#!/usr/bin/env python3
"""Audit standalone replication34 publication without importing producer code.

The four scientific payloads and archived original summary remain immutable.
A separate pinned checkout supplies historical hash inputs; they are read and
hashed, never represented as present in the publication worktree. Only current
replication34 code hashes and dependency versions may differ in the new summary.
No network, Monte Carlo simulation, participants, model fitting or CI API calls.
"""
from __future__ import annotations
import argparse
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path, PurePosixPath
import re
import subprocess
from typing import Any

HISTORY_COMMIT = '27ea33fc6cfd594c309b63ae59e695129511a9c8'
ORIGINAL_COMMIT = '661bcaf94e59a8c6a18bda0cc9771c7515b83989'
ORIGINAL_AUDITOR_SHA = '673d66c59bd35391d14135176bdd11ba26da0af518b1df4807db6e5459c9bf3d'
LOCK_SHA = '580114c2369a5277eec08dd02941c5defe2c63bb25436b2e6285389640640c23'
SOURCE_LOCK_SHA = 'd5359ff17b306e94937749915f553a7999a17244bd44314cc383efd31f4ad676'
ORIGINAL_SHA = {
    'RESULTS_SUMMARY.json': 'bbb06d8a22f7075a5e2e8fe303533f9991719510d71c13c4e4b4c93c02bb39ee',
    'PRECISION_RESULTS.json': 'ef3ab7eaf6d8a0c6fd0d7bb6232cd3241076a244a7671319e42a9c60a6e5271f',
    'SIMULATION_VALIDATION.json': 'ec9e6b0024fbb14430ec627e3c3ca55e431cfe7d286b23787652a31e3fcd2235',
    'READINESS_RESULTS.json': 'f781069d25e516ae38a7cd01abb940f6fc59eef6dcf891491bf39bef0945ca97',
    'SOURCE_MANIFEST.json': '08658be210a97b52c642f1b5277fb94af9ac4f5e723efe3ffc961668e2f35b50',
}
CURRENT_PATHS = {f'replication34/{x}' for x in ('__init__.py', 'analyze.py', 'gate.py', 'precision.py')}
PACKAGE_INITIALIZER_SHA = {
    'research/__init__.py': '83c825f1ceac051d7638bfb1b18d35df00ccf39a0125ca7d8af1811da7e4c648',
    'research/mortality/__init__.py': '4af956dca1d4405de46da9d807b39d9f9014b103e4b1fbf1f06deed56a8b18cc',
}
SUMMARY_ALLOWED = {'scientific_python_sha256', 'versions'}


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_json(path: Path) -> Any:
    # Reject duplicate keys and nonfinite JSON constants, rather than silently
    # normalizing an altered archive into the expected mathematical values.
    def unique(pairs):
        out = {}
        for key, value in pairs:
            if key in out:
                raise ValueError(f'Duplicate JSON key: {key}')
            out[key] = value
        return out
    def invalid(value):
        raise ValueError(f'Nonfinite JSON constant: {value}')
    return json.loads(path.read_text(encoding='utf-8'), object_pairs_hook=unique,
                      parse_constant=invalid)


def checked_path(root: Path, name: str) -> Path:
    relative = PurePosixPath(name)
    if relative.is_absolute() or '..' in relative.parts or not relative.parts:
        raise ValueError(f'Unsafe locked path: {name}')
    path = root.joinpath(*relative.parts)
    if not path.is_file() or path.is_symlink() or not path.resolve().is_relative_to(root.resolve()):
        raise ValueError(f'Missing, linked or escaping locked file: {name}')
    return path


def source_path(local: Path, history_mortality: Path, name: str) -> Path:
    if name.startswith('../'):
        return checked_path(history_mortality, name[3:])
    return checked_path(local, name)


def exact_tree_differences(left: Any, right: Any, prefix: str = '') -> list[str]:
    """Strict recursive compare: type, keys, order and values all matter."""
    if type(left) is not type(right):
        return [prefix or '<root>']
    if isinstance(left, dict):
        differences = [f'{prefix}/{key}' for key in left.keys() ^ right.keys()]
        for key in left.keys() & right.keys():
            differences.extend(exact_tree_differences(left[key], right[key], f'{prefix}/{key}'))
        return sorted(differences)
    if isinstance(left, list):
        if len(left) != len(right):
            return [prefix + '/length']
        return [p for index, (a, b) in enumerate(zip(left, right))
                for p in exact_tree_differences(a, b, f'{prefix}/{index}')]
    return [] if left == right else [prefix or '<root>']


def summary_differences(original: dict, current: dict) -> list[str]:
    return exact_tree_differences(
        {k: v for k, v in original.items() if k not in SUMMARY_ALLOWED},
        {k: v for k, v in current.items() if k not in SUMMARY_ALLOWED})


def load_independent_math(repo: Path):
    path = checked_path(repo, 'tools/verify_replication34_aggregates.py')
    if sha(path) != ORIGINAL_AUDITOR_SHA:
        raise ValueError('Original independent aggregate auditor changed')
    spec = importlib.util.spec_from_file_location('rep34_original_independent_math', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def git_head(root: Path) -> str:
    top = subprocess.run(['git', '-C', str(root), 'rev-parse', '--show-toplevel'],
                         check=True, capture_output=True, text=True).stdout.strip()
    if Path(top).resolve() != root.resolve():
        raise ValueError('History root is not the root of its own git checkout')
    return subprocess.run(['git', '-C', str(root), 'rev-parse', 'HEAD'], check=True,
                          capture_output=True, text=True).stdout.strip()


class Audit:
    def __init__(self):
        self.checks = []

    def check(self, name: str, passed: bool, detail: Any = None):
        row = {'name': name, 'passed': bool(passed)}
        if detail is not None:
            row['detail'] = detail
        self.checks.append(row)


def audit_math(a: Audit, independent, local: Path, data: dict):
    scenarios = read_json(local / 'SCENARIOS.json')
    simulation_cases = read_json(local / 'SIMULATION_CASES.json')
    plans = data['PRECISION_RESULTS.json']
    simulations = data['SIMULATION_VALIDATION.json']
    summary = data['RESULTS_SUMMARY.json']
    a.check('12 immutable stipulated scenarios, never empirical effect defaults',
            len(scenarios['scenarios']) == len(plans['cases']) == 12
            and scenarios['not_estimated_from_33B_effects_or_SEs'] is True
            and plans['adequate_sample_size_or_power_established'] is False)
    outputs = {}
    for given, row in zip(scenarios['scenarios'], plans['cases']):
        label = given['id']; result = row['result']
        expected = independent.expected_plan(given)
        outputs[label] = result
        a.check(label + ' exact frozen input', row['id'] == label
                and not exact_tree_differences(given, row['assumptions']))
        numeric = {k: v for k, v in expected.items() if not k.startswith('event_')}
        a.check(label + ' independent variance, normal SE, ceiling and joint yield',
                independent.subset_equal(numeric, result))
        event = result['event_count']
        a.check(label + ' raw-n binomial recurrence and event variance',
                independent.close(expected['event_expected'], event['expected_raw_events'])
                and independent.close(expected['event_variance'], event['variance_proxy'])
                and event['n_raw_trials'] == given['n_evaluable']
                and event['predictive_range'] == expected['event_quantiles']
                and event['n_over_DEFF_used_as_binomial_trials'] is False)
        floor = expected['event_probability_below_floor']
        observed = event['probability_below_technical_floor']
        a.check(label + ' exact-binomial technical-floor probability',
                floor is None and observed is None or floor is not None and observed is not None
                and math.isclose(floor, observed, rel_tol=1e-9, abs_tol=1e-300))
        a.check(label + ' no power, fitted-algorithm, censoring or recruitment overclaim',
                result['mean_effect_used_to_choose_n'] is False
                and result['algorithm_precision_status'] == 'adequacy_unqualified'
                and result['algorithm_training_uncertainty_included'] is False
                and result['full_CIF_validation_eligible'] is False
                and result['clinical_adequacy_or_power_established'] is False
                and result['IPCW_method_extension_required'] == (given['censoring_fraction'] > 0)
                and result['joint_evaluable_yield_probability'] == given['q_joint'] * given['q_known_conditional'])
    a.check('benefit/null/harm never select precision N', len({outputs[k]['conditional_n_for_half_width']
            for k in ('CENTRAL_TECHNICAL_NULL', 'HYPOTHETICAL_BENEFIT', 'HYPOTHETICAL_HARM')}) == 1)
    a.check('fourfold evaluable n halves conditional SE', independent.close(
            outputs['LARGER_EVALUABLE_DOMAIN']['conditional_SE'], outputs['CENTRAL_TECHNICAL_NULL']['conditional_SE'] / 2))
    a.check('four immutable toy laws, no production Monte Carlo rerun',
            len(simulation_cases['cases']) == len(simulations['cases']) == 4)
    moment_flags = []
    for given, row in zip(simulation_cases['cases'], simulations['cases']):
        label = given['id']; result = row['result']; expected = independent.expected_simulation(given)
        reps = given['replicates']
        a.check(label + ' exact atoms, clusters, fixed seed and 20k budget',
                row['id'] == label and not exact_tree_differences(given, row['assumptions'])
                and result['seed'] == given['seed'] and result['replicates'] == reps == 20000)
        a.check(label + ' exact rational atom moments and total variance identity',
                independent.subset_equal({k: v for k, v in expected['exact_moments'].items()
                                          if k != 'fourth_central_moment_d'}, result['exact_moments'])
                and independent.close(result['exact_moments']['variance_d'],
                    result['exact_moments']['mean_conditional_Y_variance']
                    + result['exact_moments']['variance_of_conditional_means']))
        keys = ('n_participants', 'deff', 'exact_mean_d', 'exact_variance_mean_d',
                'exact_normal_half_width', 'MCSE_mean_d', 'exact_expected_event_count',
                'exact_event_count_variance_in_this_toy', 'MCSE_mean_event_count')
        a.check(label + ' independent cluster mixture ICC, DEFF and event moments',
                independent.subset_equal({k: expected[k] for k in keys}, result))
        variance = result['MC_variance_mean_d']; fourth = result['MC_centered_fourth_moment_dbar']
        mcse = math.sqrt(max(0., (fourth - (reps - 3) / (reps - 1) * variance ** 2) / reps))
        a.check(label + ' reported fourth-moment variance MCSE', variance >= 0 and fourth >= 0
                and independent.close(result['estimated_MCSE_mean_d'], math.sqrt(variance / reps))
                and independent.close(result['approximate_fourth_moment_MCSE_variance'], mcse))
        for key in ('known_toy_variance_normal_CI_coverage', 'estimated_cluster_mean_SE_normal_CI_coverage'):
            coverage = result[key]
            a.check(label + ' ' + key + ' Wilson and binomial MCSE', coverage is not None
                    and type(coverage['covered']) is int and 0 <= coverage['covered'] <= coverage['replicates'] == reps
                    and independent.subset_equal(independent.coverage_summary(coverage['covered'], reps), coverage))
        specifications = {
            'mean_D': (result['MC_mean_d'], result['exact_mean_d'], result['MCSE_mean_d']),
            'variance_Dbar': (variance, result['exact_variance_mean_d'], mcse),
            'mean_event_count': (result['MC_mean_event_count'], result['exact_expected_event_count'], result['MCSE_mean_event_count'])}
        for key, (observed, target, standard_error) in specifications.items():
            standardized = None if standard_error == 0 else abs(observed - target) / standard_error
            passed = observed == target if standard_error == 0 else standardized <= 6
            moment_flags.append(passed); entry = row['moment_checks'][key]
            a.check(label + ' ' + key + ' predeclared six-MCSE arithmetic',
                    independent.close(entry['observed'], observed) and independent.close(entry['expected'], target)
                    and independent.close(entry['MCSE'], standard_error)
                    and independent.close(entry['absolute_discrepancy_in_MCSE'], standardized)
                    and entry['predeclared_six_MCSE_numerical_check_passed'] == passed)
        a.check(label + ' coverage reported, never forced to 95 percent',
                row['coverage_95pct_is_reported_not_required_or_tuned'] is True
                and result['positive_variance_CI_coverage_uses_no_added_tolerance'] is True
                and result['predictions_fixed_before_Y'] is True
                and result['refitted_algorithm_precision_validated'] is False
                and result['external_cohort_or_clinical_validation_performed'] is False
                and result['participant_simulation_arrays_exported'] is False)
    a.check('all 12 moment diagnostics match unchanged summary',
            all(moment_flags) == simulations['all_predeclared_moment_checks_passed'] == summary['moment_checks_passed'])


def audit_sources(a: Audit, independent, local: Path, data: dict):
    candidates = read_json(local / 'CANDIDATES.json')
    contract = read_json(local / 'GATE_CONTRACT.json')
    readiness = data['READINESS_RESULTS.json']
    sources = data['SOURCE_MANIFEST.json']
    refs = {row['id'] for row in sources['sources']}
    a.check('18 candidates and 32 unique source references',
            len(candidates) == len({c['candidate_id'] for c in candidates}) == 18
            and len(refs) == len(sources['sources']) == 32)
    required = contract['requirements']; transport = contract['transport_additional_requirements']
    a.check('306 tri-state source cells cited, no numeric laboratory/data availability inference',
            len(required) == 15 and len(transport) == 2 and all(
                set(c['checks']) == set(required + transport)
                and all(e['state'] in ('verified_true', 'verified_false', 'not_verified')
                        and all(s in refs for s in e['source_ids'])
                        and (e['state'] != 'verified_true' or bool(e['source_ids']))
                        for e in c['checks'].values()) for c in candidates))
    for mode in ('algorithm_replication', 'locked_source_transport'):
        output = readiness[mode]; requirements = required + (transport if mode == 'locked_source_transport' else [])
        a.check(mode + ' complete ordered route inventory',
                len(output['qualifications']) == output['candidate_count'] == 18
                and [x['candidate_id'] for x in output['qualifications']] == [c['candidate_id'] for c in candidates])
        ready = []
        for candidate, row in zip(candidates, output['qualifications']):
            blockers = independent.expected_qualification(candidate, requirements)
            if not blockers:
                ready.append(candidate['candidate_id'])
            a.check(mode + ' ' + candidate['candidate_id'] + ' independent fail-closed blockers',
                    row['blockers'] == blockers and row['required_checks'] == requirements
                    and row['checks'] == candidate['checks']
                    and row['source_ready'] == row['source_ready_for_protocol'] == (not blockers)
                    and row['access_actions_performed'] is False
                    and row['person_level_download_authorized'] is False
                    and row['outcome_comparison_authorized'] is False
                    and row['precision_adequacy'] == 'not_established')
        a.check(mode + ' zero qualified routes and no activation from metadata',
                ready == output['source_ready_candidate_ids'] == [] and output['qualified_candidate_count'] == 0
                and output['outcome_comparison_authorized'] is False)


def verify(repo: Path, history_root: Path, results_dir: Path, log: Path | None = None) -> dict:
    repo = repo.resolve(); history_root = history_root.resolve(); results_dir = results_dir.resolve()
    local = repo / 'research/mortality/replication34'; historical = history_root / 'research/mortality'
    a = Audit()
    independent = load_independent_math(repo)
    actual_head = git_head(history_root)
    a.check('actual historical checkout HEAD is immutable27ea', actual_head == HISTORY_COMMIT)
    # Do not continue using any other tree even if selected file hashes happen to match.
    if actual_head != HISTORY_COMMIT:
        raise ValueError('Historical checkout must be exactly ' + HISTORY_COMMIT)
    a.check('original immutable metadata lock', sha(local / 'LOCK.json') == LOCK_SHA)
    a.check('original immutable source lock', sha(local / 'SOURCE_LOCK.json') == SOURCE_LOCK_SHA)
    if sha(local / 'LOCK.json') != LOCK_SHA or sha(local / 'SOURCE_LOCK.json') != SOURCE_LOCK_SHA:
        raise ValueError('Original source/protocol lock changed')
    lock = read_json(local / 'LOCK.json'); source_lock = read_json(local / 'SOURCE_LOCK.json')
    a.check('historical lock base and inventory103 plus22', lock['base_commit'] == HISTORY_COMMIT
            and len(lock['historical_scientific_python_sha256']) == 103 and len(lock['metadata_sha256']) == 22)
    old = lock['historical_scientific_python_sha256']
    for name, digest in old.items():
        a.check('actual historical scientific bytes ' + name, sha(checked_path(historical, name)) == digest)
    for name, digest in lock['metadata_sha256'].items():
        base = local if name in lock['local_metadata'] else historical
        a.check('actual locked metadata bytes ' + name, sha(checked_path(base, name)) == digest)
    a.check('historical source tree complete scientific inventory103',
            {str(p.relative_to(historical)) for p in historical.rglob('*.py')} == set(old))
    current = {str(p.relative_to(local.parent)): sha(p) for p in sorted(local.rglob('*.py'))}
    a.check('four current runtime code paths only', set(current) == CURRENT_PATHS and not set(current) & set(old))
    for name, digest in PACKAGE_INITIALIZER_SHA.items():
        a.check('actual main runtime package initializer ' + name,
                sha(checked_path(repo, name)) == digest)
    expected_union = {**old, **current}
    original = {}; data = {}
    for name, digest in ORIGINAL_SHA.items():
        archived = checked_path(local, name)
        a.check('immutable original saved payload ' + name, sha(archived) == digest)
        original[name] = read_json(archived)
        path = checked_path(results_dir, name); data[name] = read_json(path)
        if name != 'RESULTS_SUMMARY.json':
            a.check('byte-identical reproduced scientific payload ' + name,
                    path.read_bytes() == archived.read_bytes())
    original_summary = original['RESULTS_SUMMARY.json']; summary = data['RESULTS_SUMMARY.json']
    differences = summary_differences(original_summary, summary)
    a.check('new summary permits ONLY scientific code hashes and versions', not differences, differences)
    a.check('all107 actual historical/current code hashes in publication summary',
            summary.get('scientific_python_sha256') == expected_union and len(expected_union) == 107)
    a.check('original summary retains immutable103historical plus4original hashes',
            len(original_summary['scientific_python_sha256']) == 107
            and {k: v for k, v in original_summary['scientific_python_sha256'].items() if k not in CURRENT_PATHS} == old)
    a.check('versions metadata has only explicit dependency version strings',
            set(summary.get('versions', {})) == {'python', 'numpy', 'scipy'}
            and all(isinstance(v, str) and bool(re.fullmatch(r'[0-9][A-Za-z0-9.+_-]*', v)) for v in summary['versions'].values()))
    for name, digest in source_lock['files_sha256'].items():
        a.check('reviewed local source pin ' + name, sha(checked_path(local, name)) == digest)
    sources = data['SOURCE_MANIFEST.json']
    for name, digest in sources['review_sha256'].items():
        a.check('actual local/historical source-review bytes ' + name,
                sha(source_path(local, historical, name)) == digest)
    catalog = sources['prior_catalog']
    a.check('actual pinned prior catalog bytes', sha(source_path(local, historical, catalog['path'])) == catalog['sha256'])
    a.check('old33B neutral result actual byte hash, no copied model publication',
            sha(checked_path(historical, 'cause33b/RESULTS_SUMMARY.json')) == '952a23344a6d16a453654de50bc4af437da605df9f8b2d120ced58fbbc9ec855')
    for name, digest in summary['outputs_sha256'].items():
        a.check('new summary output checksum ' + name, name in ORIGINAL_SHA and name != 'RESULTS_SUMMARY.json'
                and sha(checked_path(results_dir, name)) == digest)
    if log is not None:
        text = independent.log_text(log)
        for name in ORIGINAL_SHA:
            a.check('exact single complete stdout block ' + name,
                    independent.marked(text, name).encode('utf-8') == (results_dir / name).read_bytes())
    audit_math(a, independent, local, data)
    audit_sources(a, independent, local, data)
    failed = [x['name'] for x in a.checks if not x['passed']]
    return {
        'schema': 'replication34-publication-audit/1',
        'all_passed': not failed, 'passed': len(a.checks) - len(failed), 'total': len(a.checks),
        'original_archive_commit': ORIGINAL_COMMIT, 'historical_checkout': actual_head,
        'results_execution_origin': 'caller_supplied_files; completed_plaintext_log_checked' if log else 'caller_supplied_files; stdout_not_supplied',
        'original_aggregate_sha256': ORIGINAL_SHA,
        'reproduced_aggregate_sha256': {name: sha(results_dir / name) for name in ORIGINAL_SHA},
        'scientific_inventory': {'actual_historical_files': 103, 'actual_current_replication34_files': 4,
                                 'historical_models_imported_or_executed': False},
        'allowed_summary_changes': sorted(SUMMARY_ALLOWED),
        'independent_math_auditor_sha256': ORIGINAL_AUDITOR_SHA,
        'completed_stdout_sha256': sha(log) if log else None,
        'checks': a.checks, 'failed': failed,
        'limits': [
            'Only the four immutable scientific payloads are required byte-identical; the new summary may change only current code hashes and dependency versions.',
            'All103 historical scientific files and12 historical metadata files are read and hashed from a separate actual checkout at27ea. They are not included in or imported by the standalone iteration34 runtime.',
            'Mathematical functions are imported only from the unchanged original independent auditor, never producer code. Scenarios, rational prediction atoms, cluster DEFF, Wilson intervals and MCSE identities are recomputed; Monte Carlo simulation is not rerun by this auditor.',
            'No participant files, outcome models, external source downloads, network calls, CI API access or new cohort validation are performed. Source metadata qualification remains separate from statistical adequacy.',
            'Dependency version strings are permitted provenance metadata and checked for explicit keys and syntax; this offline auditor does not authenticate installed package versions in a different execution environment.',
            'An optional caller-supplied completed plaintext log verifies all five exported blocks; file verification alone does not claim a CI status or runtime checkout.'
        ]}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo', type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument('--history-root', type=Path, default=os.environ.get('REP34_HISTORY_ROOT'))
    parser.add_argument('--results-dir', type=Path, required=True)
    parser.add_argument('--log', type=Path)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args(argv)
    if args.history_root is None:
        parser.error('--history-root or REP34_HISTORY_ROOT is required; historical hashes cannot be assumed')
    try:
        report = verify(args.repo, args.history_root, args.results_dir, args.log)
    except (OSError, ValueError, KeyError, TypeError, subprocess.CalledProcessError) as error:
        report = {'schema': 'replication34-publication-audit/1', 'all_passed': False,
                  'passed': 0, 'total': 1, 'failed': ['audit failed closed'],
                  'error': str(error), 'checks': [{'name': 'audit failed closed', 'passed': False}]}
    args.out.write_text(json.dumps(report, ensure_ascii=False, indent=2, allow_nan=False) + '\n', encoding='utf-8')
    print(json.dumps({k: report[k] for k in ('passed', 'total', 'all_passed', 'failed')}, ensure_ascii=False))
    return 0 if report['all_passed'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
