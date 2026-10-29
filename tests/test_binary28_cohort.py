"""Synthetic clinical-definition and source-boundary tests for new binary model."""
import json,tempfile,unittest
from pathlib import Path
from unittest.mock import patch
import numpy as np
import pandas as pd
from research.mortality.binary28 import cohort,analyze

class CoreTests(unittest.TestCase):
    def frame(self):
        return pd.DataFrame(dict(SEQN=range(6),HSAGEIR=[60]*6,HSSEX=[1,2]*3,
          HAR1=[2,1,1,8,1,2],HAR3=[np.nan,1,2,1,8,np.nan],
          HAZMNK1R=[130,140,888,120,140,69],HAZNOK1R=[6,1,6,88,7,1]))
    def test_original_smoking_skip_and_unknowns(self):
        d=cohort.engineer(self.frame())
        np.testing.assert_equal(d.current[:3].to_numpy(),[0,1,0]);np.testing.assert_equal(d.former[:3].to_numpy(),[0,0,1])
        self.assertTrue(d.current.iloc[3:5].isna().all());self.assertEqual(d.current.iloc[5],0)
    def test_sbp_missing_codes_counts_and_no_frequency_clipping(self):
        d=cohort.engineer(self.frame())
        np.testing.assert_allclose(d.sbp[:2],[0,.5]);self.assertTrue(d.sbp.iloc[2:5].isna().all())
        self.assertAlmostEqual(d.sbp.iloc[5],(69-130)/20)
        f=self.frame();f.loc[0,'HAZMNK1R']=60
        self.assertAlmostEqual(cohort.engineer(f).sbp.iloc[0],(60-130)/20)
    def test_source_values_not_mutated(self):
        d=self.frame();copy=d.copy();cohort.engineer(d);pd.testing.assert_frame_equal(d,copy)
    def test_core_reader_uses_late_fields_and_two_readers(self):
        fields=json.loads((cohort.ROOT/'CLINICAL_FIELDS.json').read_text())['positions_1_based_inclusive']
        sas='LENGTH\n SEQN 7\n;\nINFILE ADULT MISSOVER LRECL=3348;\nINPUT\n'
        for name,(start,end) in fields.items():sas+=f' {name} {start}'+(f'-{end}' if end!=start else '')+'\n'
        sas+=';'
        values=dict(SEQN=1,HSSEX=2,HSAGEIR=55,HSAGEU=2,SDPPSU6=1,SDPSTRA6=12,HAR1=1,HAR3=2,HAZMNK1R=143,HAZNOK1R=6)
        chars=[' ']*3345
        for name,(start,end) in fields.items():chars[start-1:end]=str(values[name]).rjust(end-start+1)
        with tempfile.TemporaryDirectory() as t:
            root=Path(t);(root/'a.sas').write_text(sas);(root/'a.dat').write_text(''.join(chars)+'\n')
            d,a=cohort.read_core_checked(root/'a.dat',root/'a.sas')
            self.assertEqual(d.HAR3.iloc[0],2);self.assertEqual(d.HAZMNK1R.iloc[0],143)
            self.assertEqual(a['selected_fields_end'],3341);self.assertEqual(a['cell_comparisons'],10)
            (root/'a.dat').write_text(''.join(chars[:3340]))
            with self.assertRaisesRegex(ValueError,'truncation'):cohort.read_core_checked(root/'a.dat',root/'a.sas')
    def test_study_lock_works_before_any_source_parse(self):
        with patch.object(cohort,'check_study',side_effect=ValueError('locked study')):
            with patch.object(cohort.prior,'load') as source:
                with self.assertRaisesRegex(ValueError,'locked study'):cohort.load('/never-open')
                source.assert_not_called()

class DesignTests(unittest.TestCase):
    def test_partition_is_outcome_blind_and_cluster_stable(self):
        clusters=np.repeat([str(i) for i in range(20)],3)
        f=analyze.assign_outer(clusters)
        np.testing.assert_equal(analyze.assign_outer(clusters[::-1]),f[::-1])
        self.assertEqual(set(f),set(range(5)))
        for c in set(clusters):self.assertEqual(len(set(f[clusters==c])),1)
    def test_all_failed_bootstrap_has_explicit_primary_counts(self):
        r=analyze.intervals([{'status':'undefined'}]*3,3)
        for key in ('delta_brier','delta_auc'):
            self.assertIsNone(r[key]['lower']);self.assertEqual(r[key]['undefined'],3)
            self.assertTrue(r[key]['defined_draws_only'])
    def test_existing_scientific_result_not_overwritten(self):
        with tempfile.TemporaryDirectory() as t:
            with self.assertRaisesRegex(ValueError,'New output'):analyze.run('/never-open',Path(t))

if __name__=='__main__':unittest.main()
