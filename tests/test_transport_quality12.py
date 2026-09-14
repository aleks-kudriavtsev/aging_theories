"""Reporting guards and optional tests on actual new evaluation data."""
import json,os,unittest
from pathlib import Path
import numpy as np
from research.mortality.audit_transport12 import interval_status
from research.mortality.transport_no_crp12 import load_evaluation,labels,ROUTINE,sha
from research.mortality.sensitivity_assay12 import FEATURES
from research.mortality.replay_transport12 import MODEL_SHA,MODEL_PATH

class ReportingTests(unittest.TestCase):
    def test_complete_interval(self):self.assertEqual(interval_status(1000,0),'conditional_on_fixed_models')
    def test_undefined_not_nominal(self):self.assertEqual(interval_status(507,493),'defined_draws_only_not_nominal_interval')
    def test_all_undefined(self):self.assertIn('not_nominal',interval_status(0,1000))
    def test_invalid_accounting(self):
        for v,u in [(999,0),(-1,1001),(True,999)]:
            with self.assertRaises(ValueError):interval_status(v,u)
    def test_sensitivity_only_two_removed(self):self.assertEqual(set(ROUTINE)-set(FEATURES),{'rdw','log_wbc'})
    def test_primary_immutable(self):self.assertEqual(sha(MODEL_PATH),MODEL_SHA)

class DataTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        path=os.environ.get('TRANSPORT12_SOURCE')
        if not path:raise unittest.SkipTest('Actual separately retrieved sources required')
        cls.path=Path(path);cls.frames={y:load_evaluation(cls.path,y,MODEL_SHA)[0] for y in (2011,2013)}
    def test_all_source_hashes(self):
        doc=json.loads((self.path/'manifest.json').read_text())
        self.assertEqual(len(doc['files']),54)
        for r in doc['files']:self.assertEqual(sha(self.path/r['file']),r['sha256'])
    def test_count_and_disjointness(self):
        a,b=self.frames[2011],self.frames[2013];self.assertEqual((len(a),len(b)),(2666,3013));self.assertFalse(set(a.SEQN)&set(b.SEQN))
    def test_complete_domain(self):
        for d in self.frames.values():self.assertTrue(np.isfinite(d[ROUTINE]).all().all());self.assertNotIn('log_crp',d)
    def test_partition_support(self):
        for year,h,deaths in [(2011,5,130),(2013,4,113)]:
            l=labels(self.frames[year],h);self.assertEqual(int((l>0).sum()),deaths);self.assertTrue(np.isin(l,range(6)).all())

if __name__=='__main__':unittest.main()
