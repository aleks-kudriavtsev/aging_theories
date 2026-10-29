import tempfile
import unittest
from pathlib import Path
import numpy as np
import pandas as pd
from research.mortality.panel25 import analyze as a

class Panel25Tests(unittest.TestCase):
    def domain(self):
        return pd.DataFrame(dict(elig=[1]*3, dead=[0,1,0], t_exam=[60,0,60], RIAGENDR=[1,2,1], WTMEC2YR=[2.,3.,0.], WTSSNH2Y=[0.,0.,5.], SDMVSTRA=[1,1,2], SDMVPSU=[1,2,1]))
    def test_mec_not_nfl_weights(self):
        d=a.prepare_domain(self.domain());self.assertEqual(len(d),2)
        np.testing.assert_equal(d.weight.to_numpy(),[2.,3.]);self.assertEqual(d.time_years.iloc[1],.5/12)
    def test_censoring_fails_before_metrics(self):
        d=self.domain();d.loc[0,'t_exam']=47
        with self.assertRaisesRegex(ValueError,'censoring'):a.prepare_domain(d)
    def test_missing_design_fails(self):
        d=self.domain();d.loc[0,'SDMVPSU']=np.nan
        with self.assertRaises(ValueError):a.prepare_domain(d)
    def test_duplicate_or_foreign_component_ids(self):
        left=pd.DataFrame({'SEQN':[1,2]})
        for ids in ([1,1],[1,3],[1,np.nan]):
            with self.assertRaises(ValueError):a.join_component(left,pd.DataFrame({'SEQN':ids}),{1,2})
    def test_left_join_preserves_missing_profiles(self):
        d=a.join_component(pd.DataFrame({'SEQN':[1,2]}),pd.DataFrame({'SEQN':[1],'urate':[5.]}),{1,2})
        self.assertEqual(len(d),2);self.assertTrue(pd.isna(d.urate.iloc[1]))
    def arrays(self):
        rng=np.random.default_rng(25);n=100
        return np.arange(n)%3==0,np.full(n,.2),np.full(n,.25),np.ones(n),rng.normal(size=n),np.arange(n)%5
    def test_test_fold_outcomes_cannot_change_own_predictions(self):
        y,hs,h,w,x,f=self.arrays();p=a.urate_crossfit(y,hs,h,w,x,f)
        other=y.copy();other[f==0]=~other[f==0];q=a.urate_crossfit(other,hs,h,w,x,f)
        for name in p:np.testing.assert_allclose(p[name][f==0],q[name][f==0],atol=0,rtol=0)
    def test_constant_marker_matches_calibration(self):
        y,hs,h,w,x,f=self.arrays();p=a.urate_crossfit(y,hs,h,w,np.ones(len(x)),f)
        np.testing.assert_allclose(p['urate'],p['calibration'],atol=1e-12)
    def test_probabilities_and_brier(self):
        y,hs,h,w,x,f=self.arrays();p=a.urate_crossfit(y,hs,h,w,x,f)
        for v in p.values():self.assertTrue(((v>=0)&(v<=1)).all())
        s=a.score(y,p,w);expected=np.mean((y-p['urate'])**2)-np.mean((y-p['calibration'])**2)
        self.assertAlmostEqual(s['delta_brier'],expected)
    def test_undefined_bootstrap_not_silently_dropped(self):
        r=a.summarize([{'delta_brier':.01,'delta_auc':None},{'status':'undefined'}],2)
        self.assertEqual(r['delta_brier']['undefined'],1);self.assertTrue(r['delta_brier']['defined_draws_only'])
        self.assertIsNone(r['delta_auc']['lower']);self.assertEqual(r['delta_auc']['undefined'],2)
    def test_missing_sources_no_results(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);out=root/'result';(root/'manifest.json').write_text('{"files":[]}')
            with self.assertRaisesRegex(ValueError,'Missing required input'):a.run(root,out)
            self.assertFalse(out.exists())
    def test_existing_output_not_overwritten(self):
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaisesRegex(ValueError,'New output'):a.run(Path(tmp)/'missing',Path(tmp))
    def test_cluster_partition_independent_of_row_order(self):
        c=np.array([f'{j//2}:{j%2}' for j in range(20)]*2);f=a.base.assign_folds(c)
        np.testing.assert_equal(a.base.assign_folds(c[::-1]),f[::-1])
        for key in set(c):self.assertEqual(len(set(f[c==key])),1)
if __name__=='__main__':unittest.main()
