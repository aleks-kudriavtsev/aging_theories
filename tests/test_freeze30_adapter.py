"""Synthetic-only tests of the frozen reported-value research adapter.

The reference feature equations and score calculations below are independent
of the adapter's engineering and prediction helpers. No participant data is
read, scored, logged, or exported by these tests.
"""
import copy
import hashlib
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

from research.mortality.freeze30 import adapter as a


def record(**changes):
    values = dict(
        age_years=50, sex='male', smoking='current', sbp_mmhg=150.,
        sbp_measurement_count=3, bmi_kg_m2=32., bmi_weight_flag=0,
        bmi_height_flag=0, hba1c_percent=6.5, hba1c_method=1,
        creatinine_historical_mg_dl=(1. + .184) / .960,
        urine_creatinine_mg_dl=100., urine_albumin_ug_ml=10.,
        cystatin_c_reported_mg_l=.5,
        beta2_microglobulin_reported_mg_l=2.,
        beta_trace_protein_reported_mg_l=4.)
    values.update(changes)
    return values


def document(rows=None):
    contract = a.check_study()
    return dict(schema=contract['input_schema'],
                contract_sha256=a.digest(a.ROOT / 'CONTRACT.json'),
                mapping_id=contract['mapping_id'],
                units=copy.deepcopy(contract['units']),
                records=[record()] if rows is None else rows)


def synthetic_artifact():
    """Known coefficients/scaling, with no training data or fitted outcomes."""
    contract = a.check_study()
    result = dict(schema=contract['artifact_schema'], target=copy.deepcopy(contract['target']),
                  mapping_id=contract['mapping_id'],
                  contract_sha256=a.digest(a.ROOT / 'CONTRACT.json'),
                  plan_sha256=a.digest(a.ROOT / 'FREEZE_PLAN.md'), models={})
    for arm in ('routine', 'routine_all3'):
        n = 11 if arm == 'routine' else 14
        result['models'][arm] = dict(
            preprocessing=dict(model_name=arm, core_indices=list(range(11)),
                protein_transforms=[] if n == 11 else [
                    dict(column=j, kind='identity_log2') for j in range(3)],
                mean=np.linspace(-.1, .1, n).tolist(),
                scale=np.linspace(.75, 1.5, n).tolist(), constant=[False] * n),
            feature_names=list(a.CORE_NAMES) + ([] if n == 11 else list(a.PROTEIN_NAMES)),
            fit=dict(coefficients=np.linspace(-.2, .2, n).tolist(),
                     intercept=-1.25, penalty=10.), selected_penalty=10.)
    return result


def reference_features(row):
    """Direct source equations; do not call either feature-engineering module."""
    age = (row['age_years'] - 60) / 10
    female = row['sex'] == 'female'
    serum = .960 * row['creatinine_historical_mg_dl'] - .184
    ratio = serum / (.7 if female else .9)
    egfr = (142 * min(ratio, 1) ** (-.241 if female else -.302)
            * max(ratio, 1) ** -1.2 * .9938 ** row['age_years']
            * (1.012 if female else 1))
    uacr = 100 * row['urine_albumin_ug_ml'] / row['urine_creatinine_mg_dl']
    core = [age, age * age, float(not female),
            float(row['smoking'] == 'current'), float(row['smoking'] == 'former'),
            (row['sbp_mmhg'] - 130) / 20, (row['bmi_kg_m2'] - 27) / 5,
            row['hba1c_percent'], min(egfr, 60) / 15,
            max(egfr - 60, 0) / 30, np.log2(uacr)]
    proteins = [np.log2(row[name]) for name in (
        'cystatin_c_reported_mg_l', 'beta2_microglobulin_reported_mg_l',
        'beta_trace_protein_reported_mg_l')]
    return np.asarray(core), np.asarray(proteins)


class InputFeatureTests(unittest.TestCase):
    def test_source_equations_order_and_known_numeric_fixture(self):
        male = record()
        female = record(age_years=70, sex='female', smoking='former',
                        sbp_mmhg=110., bmi_kg_m2=22., hba1c_percent=5.4,
                        hba1c_method=4, creatinine_historical_mg_dl=2.,
                        urine_creatinine_mg_dl=50., urine_albumin_ug_ml=.4,
                        cystatin_c_reported_mg_l=2.,
                        beta2_microglobulin_reported_mg_l=4.,
                        beta_trace_protein_reported_mg_l=.25)
        core, proteins = a.input_features(document([male, female]))
        self.assertEqual(core.shape, (2, 11)); self.assertEqual(proteins.shape, (2, 3))
        expected = [reference_features(row) for row in (male, female)]
        np.testing.assert_allclose(core, np.vstack([x[0] for x in expected]), rtol=1e-14, atol=1e-14)
        np.testing.assert_array_equal(proteins, [[-1., 1., 2.], [1., 2., -2.]])
        # Published creatinine-only2021 equation for man50 with calibrated Scr1.
        self.assertAlmostEqual(core[0, 9], (91.69147860983348 - 60) / 30, places=12)
        self.assertEqual(core[0, 8], 4.); self.assertEqual(core[1, 9], 0.)
        self.assertAlmostEqual(core[0, 10], np.log2(10.))
        self.assertAlmostEqual(core[1, 10], np.log2(.8))

    def test_sex_smoking_and_domain_boundaries_are_explicit(self):
        rows = [record(age_years=40, smoking='never', urine_creatinine_mg_dl=10.),
                record(age_years=75, sex='female', smoking='former')]
        core, _ = a.input_features(document(rows))
        np.testing.assert_array_equal(core[:, :5], [[-2., 4., 1., 0., 0.], [1.5, 2.25, 0., 0., 1.]])

    def test_all_document_and_record_fields_are_required_and_unknown_fields_rejected(self):
        good = document()
        for key in good:
            bad = copy.deepcopy(good); del bad[key]
            with self.subTest(document_missing=key), self.assertRaises(ValueError): a.input_features(bad)
        bad = copy.deepcopy(good); bad['participant_id'] = 'not-permitted'
        with self.assertRaises(ValueError): a.input_features(bad)
        for key in good['records'][0]:
            bad = copy.deepcopy(good); del bad['records'][0][key]
            with self.subTest(record_missing=key), self.assertRaises(ValueError): a.input_features(bad)
        bad = copy.deepcopy(good); bad['records'][0]['modern_creatinine_mg_dl'] = 1.
        with self.assertRaises(ValueError): a.input_features(bad)

    def test_unknown_schema_hash_mapping_and_modern_relative_units_rejected(self):
        for key, value in [('schema', 'freeze30-input/2'), ('contract_sha256', '0' * 64),
                           ('mapping_id', 'MODERN_IDMS_STANDARDIZED_V1'),
                           ('mapping_id', 'NHANESIII_REPEAT_EXAM')]:
            bad = document(); bad[key] = value
            with self.subTest(key=key, value=value), self.assertRaises(ValueError): a.input_features(bad)
        for field, unit in [('creatinine_historical_mg_dl', 'umol/L'),
                            ('hba1c_percent', 'mmol/mol'),
                            ('cystatin_c_reported_mg_l', 'modern-standardized-mg/L'),
                            ('beta2_microglobulin_reported_mg_l', 'RFU'),
                            ('beta_trace_protein_reported_mg_l', 'NPX'),
                            ('urine_albumin_ug_ml', 'mg/L')]:
            bad = document(); bad['units'][field] = unit
            with self.subTest(field=field, unit=unit), self.assertRaises(ValueError): a.input_features(bad)
        for operation in ('missing', 'extra'):
            bad = document()
            if operation == 'missing': del bad['units']['sex']
            else: bad['units']['extra'] = 'code'
            with self.subTest(operation=operation), self.assertRaises(ValueError): a.input_features(bad)

    def test_every_numeric_field_rejects_bool_null_string_and_nonfinite(self):
        numeric = set(record()) - {'sex', 'smoking'}
        for name in sorted(numeric):
            for value in (True, False, None, '1', float('nan'), float('inf'), -float('inf')):
                with self.subTest(name=name, value=value), self.assertRaises(ValueError):
                    a.input_features(document([record(**{name: value})]))

    def test_positive_measurements_and_released_sentinels_rejected_before_formulas(self):
        positive = ('sbp_mmhg', 'bmi_kg_m2', 'hba1c_percent',
                    'creatinine_historical_mg_dl', 'urine_creatinine_mg_dl',
                    'urine_albumin_ug_ml', 'cystatin_c_reported_mg_l',
                    'beta2_microglobulin_reported_mg_l', 'beta_trace_protein_reported_mg_l')
        for name in positive:
            for value in (0., -1.):
                with self.subTest(name=name, value=value), self.assertRaises(ValueError):
                    a.input_features(document([record(**{name: value})]))
        for name, value in [('sbp_mmhg', 888), ('bmi_kg_m2', 8888),
                            ('hba1c_percent', 8888), ('creatinine_historical_mg_dl', 8888),
                            ('urine_creatinine_mg_dl', 88888), ('urine_albumin_ug_ml', 888888)]:
            with self.subTest(sentinel=name), self.assertRaises(ValueError):
                a.input_features(document([record(**{name: value})]))

    def test_integer_ranges_measurement_flags_and_categories_rejected(self):
        invalid = dict(age_years=(39, 76, 50., 50.5), sbp_measurement_count=(0, 7, 88, 3.),
                       bmi_weight_flag=(-1, 1, 2, 8, 0.), bmi_height_flag=(-1, 4, 8, 0.),
                       hba1c_method=(0, 5, 8, 1.), sex=('Male', 1, '', None),
                       smoking=('unknown', 'Never', 1, None), urine_creatinine_mg_dl=(7.9, 9.999))
        for name, values in invalid.items():
            for value in values:
                with self.subTest(name=name, value=value), self.assertRaises(ValueError):
                    a.input_features(document([record(**{name: value})]))

    def test_official_edits_methods_and_albumin_fill_are_not_clinical_qc_gates(self):
        rows = [record(bmi_weight_flag=wf, bmi_height_flag=hf, hba1c_method=method,
                       urine_albumin_ug_ml=.4)
                for wf in (0, 3, 4) for hf in (0, 1, 2, 3) for method in (1, 2, 3, 4)]
        core, proteins = a.input_features(document(rows))
        self.assertEqual(core.shape, (48, 11)); self.assertTrue(np.isfinite(proteins).all())
        np.testing.assert_array_equal(core[:, 10], np.full(48, np.log2(.4)))

    def test_nonpositive_calibrated_serum_and_derived_overflow_fail_without_flooring(self):
        for changes in ({'creatinine_historical_mg_dl': .1},
                        {'creatinine_historical_mg_dl': .184 / .960},
                        {'urine_albumin_ug_ml': 1e308}):
            with self.subTest(changes=changes), np.errstate(over='ignore'), self.assertRaises(ValueError):
                a.input_features(document([record(**changes)]))

    def test_positive_values_not_silently_clipped_to_training_frequency_range(self):
        row = record(sbp_mmhg=300., bmi_kg_m2=90., hba1c_percent=18.,
                     creatinine_historical_mg_dl=.2, urine_creatinine_mg_dl=900.,
                     urine_albumin_ug_ml=20000., cystatin_c_reported_mg_l=.0001,
                     beta2_microglobulin_reported_mg_l=1000.)
        core, proteins = a.input_features(document([row]))
        c, p = reference_features(row)
        np.testing.assert_allclose(core[0], c, rtol=1e-14, atol=1e-14)
        np.testing.assert_array_equal(proteins[0], p)

    def test_nonempty_bounded_list_and_record_type_required(self):
        for rows in ([], (), {}, None, [None], [record()] * 10001):
            bad = document(); bad['records'] = rows
            with self.subTest(kind=type(rows), size=len(rows) if rows is not None else None), self.assertRaises(ValueError):
                a.input_features(bad)
        for bad in (None, [], 'document'):
            with self.assertRaises(ValueError): a.input_features(bad)


class ArtifactAndScoringTests(unittest.TestCase):
    def test_known_frozen_scaling_coefficients_score_in_memory_without_exports_or_mutation(self):
        doc = document([record(), record(age_years=70, sex='female', smoking='never')])
        artifact = synthetic_artifact(); before_doc = copy.deepcopy(doc); before_artifact = copy.deepcopy(artifact)
        with patch.object(Path, 'write_text', side_effect=AssertionError('unexpected export')), \
             patch.object(Path, 'write_bytes', side_effect=AssertionError('unexpected export')), \
             patch('pandas.DataFrame.to_csv', side_effect=AssertionError('unexpected export')), \
             patch('numpy.savetxt', side_effect=AssertionError('unexpected export')), \
             patch('builtins.print', side_effect=AssertionError('unexpected log')):
            probabilities = a.score(doc, artifact)
        self.assertEqual(set(probabilities), {'routine', 'routine_all3'})
        features = [reference_features(row) for row in doc['records']]
        for arm, result in probabilities.items():
            spec = artifact['models'][arm]['preprocessing']; fit = artifact['models'][arm]['fit']
            raw = np.vstack([c if arm == 'routine' else np.r_[c, p] for c, p in features])
            z = (raw - np.asarray(spec['mean'])) / np.asarray(spec['scale'])
            eta = fit['intercept'] + z @ np.asarray(fit['coefficients'])
            np.testing.assert_allclose(result, 1 / (1 + np.exp(-eta)), rtol=1e-14, atol=1e-14)
            self.assertEqual(result.shape, (2,)); self.assertTrue(((result > 0) & (result < 1)).all())
        self.assertEqual(doc, before_doc); self.assertEqual(artifact, before_artifact)

    def test_routine_arm_ignores_protein_values_but_still_requires_complete_input(self):
        artifact = synthetic_artifact()
        p = a.score(document(), artifact)
        q = a.score(document([record(beta_trace_protein_reported_mg_l=16.)]), artifact)
        np.testing.assert_array_equal(p['routine'], q['routine'])
        self.assertNotEqual(p['routine_all3'][0], q['routine_all3'][0])
        bad = document(); del bad['records'][0]['beta_trace_protein_reported_mg_l']
        with self.assertRaises(ValueError): a.score(bad, artifact)

    def test_artifact_contract_hash_plan_hash_and_exact_arm_registry(self):
        artifact = synthetic_artifact(); self.assertIs(a.validate_artifact(artifact), artifact)
        for key, value in [('schema', 'freeze30-model/2'), ('contract_sha256', '0' * 64),
                           ('plan_sha256', '0' * 64), ('models', {}), ('models', [])]:
            bad = copy.deepcopy(artifact); bad[key] = value
            with self.subTest(key=key), self.assertRaises(ValueError): a.validate_artifact(bad)
        bad = copy.deepcopy(artifact); bad['models']['rank_all3'] = copy.deepcopy(bad['models']['routine_all3'])
        with self.assertRaises(ValueError): a.validate_artifact(bad)
        for bad in (None, [], 'artifact'):
            with self.assertRaises(ValueError): a.validate_artifact(bad)

    def test_changed_feature_order_selection_or_protein_transform_is_rejected(self):
        for arm in ('routine', 'routine_all3'):
            for field, value in [('model_name', 'rank_all3'), ('core_indices', list(range(10))),
                                 ('core_indices', list(reversed(range(11))))]:
                bad = synthetic_artifact(); bad['models'][arm]['preprocessing'][field] = value
                with self.subTest(arm=arm, field=field), self.assertRaises(ValueError): a.validate_artifact(bad)
            bad = synthetic_artifact(); bad['models'][arm]['feature_names'].reverse()
            with self.subTest(arm=arm), self.assertRaises(ValueError): a.validate_artifact(bad)
        for transforms in ([], [{'column': 0, 'kind': 'identity_log2'}],
                           [{'column': j, 'kind': 'weighted_midpoint_tail_clip', 'lower': -1., 'upper': 1.} for j in range(3)],
                           [{'column': j, 'kind': 'identity_log2', 'unapproved': True} for j in range(3)],
                           [{'column': j, 'kind': 'identity_log2'} for j in (2, 1, 0)]):
            bad = synthetic_artifact(); bad['models']['routine_all3']['preprocessing']['protein_transforms'] = transforms
            with self.subTest(transforms=transforms), self.assertRaises(ValueError): a.validate_artifact(bad)
        bad = synthetic_artifact(); bad['models']['routine']['preprocessing']['protein_transforms'] = [dict(column=0, kind='identity_log2')]
        with self.assertRaises(ValueError): a.validate_artifact(bad)

    def test_all_frozen_numeric_vectors_require_exact_shape_finite_numbers_and_positive_scale(self):
        for arm in ('routine', 'routine_all3'):
            n = 11 if arm == 'routine' else 14
            for parent, field in [('preprocessing', 'mean'), ('preprocessing', 'scale'), ('fit', 'coefficients')]:
                for value in ([0.] * (n - 1), [0.] * (n + 1), [[0.]] * n,
                              [False] * n, [None] * n, ['1'] * n,
                              [float('nan')] * n, [float('inf')] * n, tuple([1.] * n)):
                    bad = synthetic_artifact(); bad['models'][arm][parent][field] = value
                    with self.subTest(arm=arm, field=field, value_type=type(value)), self.assertRaises(ValueError):
                        a.validate_artifact(bad)
            for value in (0., -1.):
                bad = synthetic_artifact(); bad['models'][arm]['preprocessing']['scale'][0] = value
                with self.subTest(arm=arm, scale=value), self.assertRaises(ValueError): a.validate_artifact(bad)
            for value in ([False] * (n - 1), [0] * n, [None] * n):
                bad = synthetic_artifact(); bad['models'][arm]['preprocessing']['constant'] = value
                with self.subTest(arm=arm, constant=value), self.assertRaises(ValueError): a.validate_artifact(bad)

    def test_fit_intercept_penalty_and_selected_penalty_are_checked(self):
        for key, values in [('intercept', (True, None, '0', float('nan'), float('inf'))),
                            ('penalty', (True, None, '10', 0., -1., 2., float('nan'), float('inf')))]:
            for value in values:
                bad = synthetic_artifact(); bad['models']['routine']['fit'][key] = value
                with self.subTest(key=key, value=value), self.assertRaises(ValueError): a.validate_artifact(bad)
        bad = synthetic_artifact(); bad['models']['routine']['selected_penalty'] = 100.
        with self.assertRaises(ValueError): a.validate_artifact(bad)

    def test_bool_must_not_alias_frozen_feature_indices_or_selected_penalty(self):
        for location in ('core_indices', 'protein_column', 'selected_penalty'):
            bad = synthetic_artifact()
            if location == 'core_indices':
                bad['models']['routine']['preprocessing']['core_indices'][0] = False
            elif location == 'protein_column':
                bad['models']['routine_all3']['preprocessing']['protein_transforms'][0]['column'] = False
            else:
                bad['models']['routine']['fit']['penalty'] = 1.
                bad['models']['routine']['selected_penalty'] = True
            with self.subTest(location=location), self.assertRaises(ValueError): a.validate_artifact(bad)


class JsonAndIntegrityTests(unittest.TestCase):
    def test_strict_json_rejects_duplicate_keys_at_every_depth_and_nonfinite_tokens(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'metadata.json'
            for payload in ('{"x":1,"x":2}', '{"outer":{"x":1,"x":2}}',
                            '{"x":NaN}', '{"x":Infinity}', '{"x":-Infinity}'):
                path.write_text(payload)
                with self.subTest(payload=payload), self.assertRaises(ValueError): a.strict_json(path)

    def test_strict_json_rejects_numeric_exponent_overflow(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'overflow.json'
            path.write_text('{"nested":[{"x":1e999}]}')
            with self.assertRaises(ValueError): a.strict_json(path)

    def test_strict_json_size_bound_checked_before_parsing(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'large.json'
            path.write_bytes(b' ' * 2_000_001)
            with self.assertRaisesRegex(ValueError, 'size limit'): a.strict_json(path)

    def test_manifest_pins_exact_artifact_bytes_and_study_hashes(self):
        artifact = synthetic_artifact()
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'model.json'; manifest_path = Path(directory) / 'manifest.json'
            path.write_text(json.dumps(artifact, sort_keys=True, allow_nan=False))
            manifest = dict(schema='freeze30-manifest/1',
                artifact_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                contract_sha256=artifact['contract_sha256'], plan_sha256=artifact['plan_sha256'])
            manifest_path.write_text(json.dumps(manifest))
            self.assertEqual(a.load_artifact(path, manifest_path), artifact)
            for key in manifest:
                bad = copy.deepcopy(manifest); bad[key] = 'wrong'
                manifest_path.write_text(json.dumps(bad))
                with self.subTest(key=key), self.assertRaises(ValueError): a.load_artifact(path, manifest_path)
            manifest_path.write_text(json.dumps(manifest)); original = path.read_bytes()
            path.write_bytes(original + b'\n')
            with self.assertRaisesRegex(ValueError, 'manifest mismatch'): a.load_artifact(path, manifest_path)

    def test_matching_byte_hash_cannot_authorize_unapproved_transforms_or_duplicate_json(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'model.json'; manifest_path = Path(directory) / 'manifest.json'
            artifact = synthetic_artifact()
            artifact['models']['routine_all3']['preprocessing']['protein_transforms'][0]['kind'] = 'modern_standardized'
            for payload in (json.dumps(artifact), '{"schema":"freeze30-model/1","schema":"freeze30-model/1"}'):
                path.write_text(payload)
                manifest_path.write_text(json.dumps(dict(schema='freeze30-manifest/1',
                    artifact_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                    contract_sha256=a.digest(a.ROOT / 'CONTRACT.json'),
                    plan_sha256=a.digest(a.ROOT / 'FREEZE_PLAN.md'))))
                with self.subTest(payload_kind=payload[:30]), self.assertRaises(ValueError):
                    a.load_artifact(path, manifest_path)

    def test_duplicate_manifest_key_rejected_even_when_hash_is_correct(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'model.json'; manifest_path = Path(directory) / 'manifest.json'
            artifact = synthetic_artifact(); path.write_text(json.dumps(artifact))
            manifest = dict(schema='freeze30-manifest/1', artifact_sha256=a.digest(path),
                            contract_sha256=artifact['contract_sha256'], plan_sha256=artifact['plan_sha256'])
            payload = json.dumps(manifest)
            payload = payload[:-1] + ',"artifact_sha256":' + json.dumps(manifest['artifact_sha256']) + '}'
            manifest_path.write_text(payload)
            with self.assertRaisesRegex(ValueError, 'Duplicate'): a.load_artifact(path, manifest_path)

    def test_study_lock_rejects_changed_local_metadata_without_touching_repository(self):
        original_root = a.ROOT
        lock = a.strict_json(original_root / 'LOCK.json')
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory) / 'freeze30'; root.mkdir()
            (root / 'LOCK.json').write_text(json.dumps(lock))
            for name in lock['metadata_sha256']:
                source = original_root.parent / name if '/' in name else original_root / name
                target = root.parent / name if '/' in name else root / name
                target.parent.mkdir(parents=True, exist_ok=True); target.write_bytes(source.read_bytes())
            with patch.object(a, 'ROOT', root):
                self.assertEqual(a.check_study()['schema'], 'freeze30-contract/1')
                for name in lock['metadata_sha256']:
                    target = root.parent / name if '/' in name else root / name
                    original = target.read_bytes(); target.write_bytes(original + b'\n')
                    with self.subTest(name=name), self.assertRaisesRegex(ValueError, 'metadata changed'): a.check_study()
                    target.write_bytes(original)


if __name__ == '__main__':
    unittest.main()
