"""Human-operated, outcome-blind joint-profile audit; standard library only.

Input is a holder-prepared presence table, itself HRS microdata. It is NOT to be
uploaded to ChatGPT, AI IDEs, GitHub or Vercel. The program never downloads data,
fits coefficients, estimates mortality, or exports person-level records.
"""
from __future__ import annotations
import argparse
import csv
import json
import re
from pathlib import Path
from .access import check_access, strict_json

PROTEINS = {'NT_PROBNP': 'P16860', 'CST3': 'P01034', 'GDF15': 'Q99988',
            'IL6': 'P05231', 'TNF': 'P01375', 'GFAP': 'P14136', 'NEFL': 'P07196'}
CORE13 = ('total_cholesterol', 'hdl', 'hba1c', 'creatinine', 'uacr', 'albumin',
          'rdw_cv', 'wbc', 'alp', 'ggt', 'bilirubin_total', 'platelets', 'nlr')
CLINICAL = ('age', 'sex', 'smoking', 'bmi', 'sbp', 'bp_treatment',
            'diabetes_history', 'cvd_history', 'cancer_history')
ALL_FIELDS = set(PROTEINS) | set(CORE13) | set(CLINICAL)
COLUMNS = ('HHID', 'PN', 'baseline_token', 'component', 'usable')


def check_mapping(mappings, observed):
    if not isinstance(mappings, dict):
        raise ValueError('mapping_object_required')
    for k in observed:
        m = mappings.get(k)
        if not isinstance(m, dict) or m.get('reviewed') is not True:
            raise ValueError('unreviewed_component_mapping')
        if not all(isinstance(m.get(f), str) and m[f].strip()
                   for f in ('source_product', 'source_variable', 'source_locator')):
            raise ValueError('mapping_provenance_missing')
        if k in PROTEINS and m.get('uniprot') != PROTEINS[k]:
            raise ValueError('protein_identity_mismatch')
        if k == 'TNF' and 'TNFR' in m['source_variable'].upper():
            raise ValueError('receptor_not_direct_TNF')
        if k == 'NT_PROBNP' and m.get('measured_form') != 'NT-proBNP':
            raise ValueError('NT_proBNP_form_not_confirmed')
        if k in PROTEINS or k in CORE13:
            if m.get('assay_scale_and_QC_reviewed') is not True:
                raise ValueError('assay_scale_and_QC_unreviewed')


def summarize_records(records, config):
    gate = check_access(config)
    if gate['status'] != 'declared_conditions_met':
        raise ValueError('access_check_blocked')
    if config.get('data_class') not in ('synthetic_fixture', 'hrs_microdata'):
        raise ValueError('record_input_wrong_data_class')
    anchor = config.get('baseline_token')
    if not isinstance(anchor, str) or not anchor.strip():
        raise ValueError('reviewed_baseline_token_required')
    common = config.get('common_core_fields')
    if not isinstance(common, list) or not common or len(common) != len(set(common)):
        raise ValueError('explicit_unique_common_core_required')
    if any(k not in CORE13 for k in common):
        raise ValueError('common_core_must_use_known_routine_components')
    people, seen, observed = {}, set(), set()
    for r in records:
        if not isinstance(r, dict) or set(r) != set(COLUMNS):
            raise ValueError('presence_schema_mismatch')
        if not isinstance(r['HHID'], str) or re.fullmatch(r'\d{6}', r['HHID']) is None:
            raise ValueError('six_character_HHID_required')
        if not isinstance(r['PN'], str) or re.fullmatch(r'\d{3}', r['PN']) is None:
            raise ValueError('three_character_PN_required')
        if r['baseline_token'] != anchor:
            raise ValueError('different_baseline_not_joined')
        k = r['component']
        if not isinstance(k, str) or k not in ALL_FIELDS:
            raise ValueError('unknown_component')
        if r['usable'] not in ('0', '1'):
            raise ValueError('usable_requires_0_or_1')
        person = (r['HHID'], r['PN'])
        key = person + (k,)
        if key in seen:
            raise ValueError('duplicate_person_component_not_averaged')
        seen.add(key)
        observed.add(k)
        people.setdefault(person, set())
        if r['usable'] == '1':
            people[person].add(k)
    if not people:
        raise ValueError('empty_presence_table')
    check_mapping(config.get('mappings'), observed)
    clinical = set(CLINICAL)
    six = set(PROTEINS) - {'TNF'}
    core = clinical | set(common)
    groups = {'clinical': clinical, 'clinical_core13': clinical | set(CORE13),
              'clinical_common_core': core,
              'common_core_plus_GDF15_IL6_GFAP': core | {'GDF15', 'IL6', 'GFAP'},
              'common_core_plus6_without_direct_TNF': core | six,
              'common_core_plus7_including_direct_TNF': core | set(PROTEINS)}
    counts = {name: sum(required <= available for available in people.values())
              for name, required in groups.items()}
    return {'status': 'local_feasibility_only', 'data_class': config['data_class'],
            'total_participants_in_presence_table': len(people), 'counts': counts,
            'common_core_fields': common,
            'absent_components': sorted(ALL_FIELDS - observed),
            'participant_records_exported': False, 'outcomes_read': False,
            'coefficients_fitted': False, 'existing_model_applicable': False,
            'approved_for_sharing': False,
            'warning': 'INTERNAL ONLY: overlapping counts require human disclosure review; do not upload this file to AI.'}


def review_partition(counts, reviewed=False):
    """Helper for ONE mutually exclusive table. Not general disclosure control.

Suppress the whole table when any cell is below five; do not publish totals or
complements then. Other publications/differencing still require human review.
"""
    if not isinstance(counts, dict) or not counts or any(
            type(n) is not int or n < 0 for n in counts.values()):
        raise ValueError('nonnegative_integer_partition_counts_required')
    if reviewed is not True or any(n < 5 for n in counts.values()):
        return {'status': 'withheld', 'counts': None, 'total': None}
    return {'status': 'reviewed_single_partition', 'counts': dict(counts),
            'total': sum(counts.values()), 'cross_release_disclosure_verified': False}


def run(config_path, presence_path, out):
    cpath, src, dest = map(Path, (config_path, presence_path, out))
    if cpath.stat().st_size > 1_000_000:
        raise ValueError('config_size_limit')
    config = strict_json(cpath.read_text(encoding='utf-8'))
    gate = check_access(config)
    if gate['status'] != 'declared_conditions_met':
        return gate  # Deliberately no stat/open of person-level source before gate.
    if config.get('data_class') not in ('synthetic_fixture', 'hrs_microdata'):
        return {'status': 'blocked', 'blockers': ['record_input_wrong_data_class']}
    if dest.exists() or dest.is_symlink():
        raise ValueError('new_output_path_required')
    if src.stat().st_size > 20_000_000:
        raise ValueError('presence_size_limit')
    with src.open(newline='', encoding='utf-8-sig') as f:
        reader = csv.DictReader(f)
        if reader.fieldnames != list(COLUMNS):
            raise ValueError('exact_presence_header_required')
        result = summarize_records(reader, config)
    with dest.open('x', encoding='utf-8') as f:
        json.dump(result, f, ensure_ascii=False, indent=2, allow_nan=False)
    return result


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--config', type=Path, required=True)
    p.add_argument('--presence', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    a = p.parse_args()
    try:
        result = run(a.config, a.presence, a.out)
        # No exact counts or identifiers written to console logs.
        print(json.dumps({'status': result['status'], 'approved_for_sharing': False}))
        raise SystemExit(2 if result['status'] == 'blocked' else 0)
    except (ValueError, OSError, UnicodeError, KeyError, TypeError):
        print(json.dumps({'status': 'failed', 'error': 'input_or_local_environment_check_failed'}))
        raise SystemExit(2)
