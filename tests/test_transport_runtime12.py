"""Autonomous replay tests: no scientific packages or participant records."""
from copy import deepcopy
import hashlib,math,tempfile,unittest
from pathlib import Path
from research.mortality.replay_transport12 import example,run,load_bundle,MODEL_PATH,MODEL_SHA,score_engineered,GROUPS

class StandaloneTests(unittest.TestCase):
    def test_checksum(self):self.assertEqual(hashlib.sha256(MODEL_PATH.read_bytes()).hexdigest(),MODEL_SHA)
    def test_synthetic(self):self.assertEqual(run(example())['declaration'],'synthetic_not_a_person')
    def test_no_imputation(self):self.assertFalse(run(example())['imputation_applied'])
    def test_no_clinical_use(self):self.assertFalse(run(example())['clinical_use_ready'])
    def test_no_noninfectious_substitution(self):self.assertIsNone(run(example())['noninfectious_natural_probability'])
    def test_no_crp(self):self.assertNotIn('crp',example()['measurements'])
    def test_missing_values(self):
        for key in example()['measurements']:
            p=example();p['measurements'].pop(key);self.assertEqual(run(p)['status'],'blocked')
    def test_invalid_country(self):
        p=example();p['country']='RU';self.assertEqual(run(p)['status'],'blocked')
    def test_horizon(self):
        for h in ([20],[8],[True],[],[1,1]):
            p=example();p['horizons_years']=h;self.assertEqual(run(p)['status'],'blocked')
    def test_mass(self):
        for r in run(example())['probabilities']:self.assertAlmostEqual(r['survival']+sum(r['cause_specific_cif'].values()),1)
    def test_monotonicity(self):
        rows=run(example())['probabilities']
        for g in GROUPS:self.assertEqual(sorted(r['cause_specific_cif'][g] for r in rows),[r['cause_specific_cif'][g] for r in rows])
    def test_identity_not_allowed(self):
        p=example();p['patient_name']='SYNTHETIC';self.assertEqual(run(p)['status'],'blocked')
    def test_finite_handling(self):
        for v in (True,math.inf,math.nan,10**400):
            p=example();p['measurements']['creatinine']['value']=v;self.assertEqual(run(p)['status'],'blocked')
    def test_units(self):
        p=example();p['measurements']['creatinine']['unit']='wrong';self.assertEqual(run(p)['status'],'blocked')
    def test_material(self):
        p=example();p['measurements']['uacr']['matrix']='serum';self.assertEqual(run(p)['status'],'blocked')
    def test_tamper(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'bad.json';p.write_text('{}')
            with self.assertRaises(ValueError):run(example(),bundle_path=p)
    def test_unmodified_input(self):
        p=example();saved=deepcopy(p);run(p);self.assertEqual(p,saved)
    def test_clinical_panel(self):self.assertEqual(run(example('clinical'))['status'],'calculated_research_only')

if __name__=='__main__':unittest.main()
