"""Numerical and actual-source tests for amended compact-panel evaluation."""
import json,os,unittest
from pathlib import Path
import numpy as np
from research.mortality.compact_panel14 import (aggregate_states,late_labels,comparator_models,sha,ROUTINE,NAMES)
from research.mortality.evaluate_compact14 import load_new
from research.mortality.ipcw14 import outcome_weights,brier_states,binary_scores,aalen_johansen_reference
from research.mortality.replay_compact14 import MODEL_SHA,MODEL_PATH

class IPCWTests(unittest.TestCase):
    def test_no_censoring_identity(self):
        o=outcome_weights([1,0,1],[1,np.nan,2],[1,5,4],[1,2,3],3)
        np.testing.assert_array_equal(o['effective_weights'],[1,2,3]);self.assertEqual(o['early_censored_n'],0)
    def test_early_censor_is_unknown_not_negative(self):
        o=outcome_weights([1,0,1],[1,np.nan,2],[1,2,3],[1,1,1],3)
        np.testing.assert_array_equal(o['effective_weights'],[1,0,2]);self.assertFalse(o['known'][1]);self.assertEqual(o['early_censored_n'],1)
    def test_censor_at_horizon_known(self):
        o=outcome_weights([0,0],[np.nan,np.nan],[3,4],[1,1],3)
        np.testing.assert_array_equal(o['effective_weights'],[1,1])
    def test_death_at_horizon_event(self):
        o=outcome_weights([1,0],[2,np.nan],[3,4],[1,1],3);self.assertEqual(o['target'][0],2)
    def test_unknown_cause_is_event(self):
        o=outcome_weights([1,0],[np.nan,np.nan],[1,4],[1,1],3);self.assertEqual(o['target'][0],3)
    def test_late_coarse_only(self):
        with self.assertRaises(ValueError):outcome_weights([1],[5],[1],[1],3)
    def test_no_horizon_support(self):
        with self.assertRaises(ValueError):outcome_weights([0],[np.nan],[1],[1],3)
    def test_input_validation(self):
        examples=[([],[],[],[],3),([0],[np.nan],[4],[-1],3),([1],[1],[-1],[1],3),([2],[1],[4],[1],3),([0],[np.inf],[4],[1],3),([0],[1],[4],[1],3),([0],[np.nan],[np.nan],[1],3),([0],[np.nan],[4],[1],True)]
        for args in examples:
            with self.subTest(args=args),self.assertRaises(ValueError):outcome_weights(*args)
    def test_mass_with_tied_death_and_censor(self):
        d=[1,0,1,0];c=[1,np.nan,2,np.nan];t=[1,1,2,4];w=[2,3,4,5]
        o=outcome_weights(d,c,t,w,3);self.assertAlmostEqual(o['effective_weights'].sum(),sum(w))
        np.testing.assert_allclose(np.bincount(o['target'],weights=o['effective_weights'],minlength=4)/sum(w),aalen_johansen_reference(d,c,t,w,3),atol=1e-14)
    def test_repeated_random_AJ_agreement(self):
        rng=np.random.default_rng(15)
        for _ in range(100):
            d=rng.integers(0,2,50);t=rng.integers(1,8,50).astype(float);w=rng.uniform(.1,3,50)
            c=np.where(d==1,rng.choice([1,2,10],50),np.nan)
            o=outcome_weights(d,c,t,w,3)
            np.testing.assert_allclose(np.bincount(o['target'],weights=o['effective_weights'],minlength=4)/sum(w),aalen_johansen_reference(d,c,t,w,3),atol=1e-13)
    def test_weight_scaling_invariance(self):
        a=outcome_weights([1,0,1],[1,np.nan,2],[1,2,3],[1,2,3],3)
        b=outcome_weights([1,0,1],[1,np.nan,2],[1,2,3],[10,20,30],3)
        np.testing.assert_allclose(a['ipcw_factors'],b['ipcw_factors'])
    def test_four_state_brier_hand(self):
        target=np.array([0,1]);p=np.array([[.8,.2,0,0],[.4,.6,0,0]])
        self.assertAlmostEqual(brier_states(target,p,[1,1],2),(.08+.32)/2)
    def test_ipcw_full_denominator(self):
        o=outcome_weights([1,0,1],[1,np.nan,2],[1,2,3],[1,1,1],3)
        p=np.full((3,4),.25)
        self.assertAlmostEqual(brier_states(o['target'],p,o['effective_weights'],3),.75)
    def test_predictions_average_keeps_censored(self):
        o=outcome_weights([1,0,1],[1,np.nan,2],[1,2,3],[1,1,1],3)
        s=binary_scores(o['target'],[.1,.9,.2],[1,1,1],o)
        self.assertAlmostEqual(s['predicted'],.4);self.assertAlmostEqual(s['observed'],1)
    def test_undefined_auc_retained(self):
        o=outcome_weights([0,0],[np.nan,np.nan],[4,5],[1,1],3)
        self.assertIsNone(binary_scores(o['target'],[.1,.2],[1,1],o)['auc'])
    def test_zero_bootstrap_weights(self):
        o=outcome_weights([1,0,1],[1,np.nan,2],[1,2,4],[1,0,1],3)
        np.testing.assert_array_equal(o['effective_weights'],[1,0,1])
    def test_invalid_probability_mass(self):
        with self.assertRaises(ValueError):brier_states([0],[[.9,.9]],[1],1)
    def test_aggregation_does_not_lose_other_causes(self):
        p=aggregate_states([[.5,.1,.1,.1,.1,.1]])
        np.testing.assert_allclose(p,[[.5,.1,.1,.3]])
    def test_aggregation_invalid(self):
        for p in [[[1,1,1,1,1,1]],[[.5,.5]],[[1,0,0,0,0,np.nan]]]:
            with self.assertRaises(ValueError):aggregate_states(p)
    def test_original_lock_guard_still_stops(self):
        with self.assertRaisesRegex(ValueError,'Early censoring'):late_labels([0],[np.nan],[35/12],3)
    def test_original_training_code_preserved(self):
        model=json.loads(MODEL_PATH.read_text())
        p=MODEL_PATH.parents[1]/'compact_panel14.py'
        self.assertEqual(sha(p),model['training_code_sha256'])
    def test_comparator_versions_unchanged(self):self.assertEqual(set(comparator_models()),{'clinical','six6','full8'})

class ActualSources(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        source=os.environ.get('NHANES_COMPACT14_SOURCE')
        if not source:raise unittest.SkipTest('Actual2015 source archive required')
        cls.source=Path(source);cls.d,cls.flow=load_new(cls.source,MODEL_SHA)
    def test_verified_files(self):self.assertEqual(self.flow['verified_input_files'],27)
    def test_reader_cells(self):self.assertEqual(self.flow['reader_audit']['cell_comparisons'],79768)
    def test_complete_common_domain(self):self.assertEqual(int(np.isfinite(self.d[ROUTINE]).all(axis=1).sum()),2922)
    def test_censoring_explicit(self):self.assertEqual(int((self.d.dead.eq(0)&self.d.t_exam.lt(36)).sum()),4)
    def test_population_flow(self):self.assertEqual(len(self.d),3256)
    def test_source_LOD_flags_retained(self):
        from research.mortality.xpt_checked import read_xpt_checked
        a=read_xpt_checked(self.source/'2015_ALB_CR_I.xpt')
        self.assertIn('URDUMALC',a);self.assertIn('URDUCRLC',a)

if __name__=='__main__':unittest.main()
