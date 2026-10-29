"""Synthetic orchestration/source-boundary tests, no real cohort records."""
import tempfile,unittest
from pathlib import Path
from unittest.mock import patch
import numpy as np
from research.mortality.routine29 import analyze,cohort,retrieve

class DesignTests(unittest.TestCase):
    def test_partition_stress_is_cluster_stable_outcome_blind(self):
        clusters=np.repeat([str(i) for i in range(20)],3)
        for s in range(5):
            f=analyze.partition(clusters,s)
            np.testing.assert_equal(analyze.partition(clusters[::-1],s),f[::-1])
            for c in set(clusters):self.assertEqual(len(set(f[clusters==c])),1)
        self.assertFalse(np.array_equal(analyze.partition(clusters,0),analyze.partition(clusters,1)))
    def test_preflight_denial_stops_before_source_manifest(self):
        with patch.object(cohort,'check_study',side_effect=ValueError('locked study')):
            with patch.object(Path,'read_text') as read:
                with self.assertRaisesRegex(ValueError,'locked study'):cohort.verify_sources('/never-open')
                read.assert_not_called()
    def test_output_exists_prevents_any_recalculation(self):
        with tempfile.TemporaryDirectory() as t:
            with patch.object(analyze,'load') as source:
                with self.assertRaisesRegex(ValueError,'New output'):analyze.run('/never-open',Path(t))
                source.assert_not_called()
    def test_all_bootstrap_failures_still_count_primary_undefined(self):
        values=analyze.intervals([dict(status='undefined')]*2,2)
        for key in ('delta_brier','delta_auc'):
            self.assertIsNone(values[key]['lower']);self.assertEqual(values[key]['undefined'],2)
            self.assertTrue(values[key]['defined_draws_only'])

class ReleaseTests(unittest.TestCase):
    def test_resolve_original_lab_exam_and_no_repeat(self):
        h='<a href="/nchs/data/nhanes3/1a/LAB.DAT">L</a><a href="/nchs/data/nhanes3/1a/exam.dat">E</a><a href="/nchs/data/nhanes3/3a/examse.dat">repeat</a>'
        r=retrieve.resolve_additional(h);self.assertEqual(set(r),{'lab.dat','exam.dat'})
        self.assertIn('/1a/',r['exam.dat'])
    def test_repeat_or_external_lookalike_cannot_supply_exam(self):
        for url in ('/nchs/data/nhanes3/3a/examse.dat','https://wwwn.cdc.gov.evil.test/nchs/data/nhanes3/1a/exam.dat'):
            html='<a href="/nchs/data/nhanes3/1a/lab.dat">L</a><a href="'+url+'">E</a>'
            with self.assertRaisesRegex(ValueError,'Unresolved'):retrieve.resolve_additional(html)
    def test_retrieval_checks_study_before_network(self):
        with patch.object(retrieve,'check_study',side_effect=ValueError('locked')):
            with patch.object(retrieve.prior,'retrieve') as network:
                with self.assertRaisesRegex(ValueError,'locked'):retrieve.retrieve('/never-create')
                network.assert_not_called()

if __name__=='__main__':unittest.main()
