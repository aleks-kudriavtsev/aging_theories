"""Synthetic-only routine29 feature-selection/nesting tests; no cohort input."""
import copy
import unittest
from unittest.mock import patch
import numpy as np
from research.mortality.binary28 import model as old
from research.mortality.routine29 import model as m


def synthetic(seed=29, n=300):
    rng=np.random.default_rng(seed)
    core=rng.normal(size=(n,11)); proteins=rng.normal(size=(n,3))
    y=(rng.random(n)<.3).astype(int); w=rng.lognormal(0,.4,n)
    clusters=np.array([f'psu{i//5}' for i in range(n)])
    folds=m.psu_partition(clusters,'binary28:outer',5)
    return core,proteins,y,w,clusters,folds


class PreprocessingTests(unittest.TestCase):
    def test_simple_and_routine_use_exact_registered_columns(self):
        core,proteins,y,w,*_=synthetic(n=100)
        for name,nfeatures in [('simple',6),('simple_all3',9),('routine',11),
                              ('routine_all3',14),('rank_all3',14),('clip_all3',14)]:
            with self.subTest(name=name):
                fit=m.fit_preprocessing(core,proteins,w,name)
                transformed=m.transform_features(core,proteins,fit)
                self.assertEqual(transformed.shape,(100,nfeatures))
                self.assertEqual(fit['core_indices'],list(range(6 if name.startswith('simple') else 11)))

    def test_simple_excludes_extra_routine_and_unused_protein_values(self):
        core,proteins,_,w,*_=synthetic(n=100)
        changed=core.copy();changed[:,6:]+=1e8
        for name in ('simple','simple_all3'):
            a=m.fit_preprocessing(core,proteins,w,name)
            b=m.fit_preprocessing(changed,proteins,w,name)
            self.assertEqual(a,b)
            np.testing.assert_array_equal(m.transform_features(core,proteins,a),m.transform_features(changed,proteins,b))
        a=m.fit_preprocessing(core,proteins,w,'simple')
        b=m.fit_preprocessing(core,proteins+1000,w,'simple')
        self.assertEqual(a,b)

    def test_invalid_shapes_missing_values_and_unregistered_arms_fail(self):
        core,proteins,_,w,*_=synthetic(n=100)
        for c,p,name in [(core[:,:6],proteins,'routine'),(core,proteins[:,:2],'routine'),
                         (core,proteins,'core')]:
            with self.assertRaises(ValueError):m.fit_preprocessing(c,p,w,name)
        core[0,8]=np.nan
        with self.assertRaises(ValueError):m.fit_preprocessing(core,proteins,w,'simple')

    def test_constant_and_zero_weight_routine_columns_use_training_support(self):
        core=np.zeros((4,11));core[:,8]=[1,2,2,1e8]
        proteins=np.ones((4,3));w=np.array([1.,2.,1.,0.])
        fit=m.fit_preprocessing(core,proteins,w,'routine_all3')
        self.assertAlmostEqual(fit['mean'][8],1.75)
        self.assertFalse(fit['constant'][8]);self.assertTrue(fit['constant'][0])
        self.assertEqual(fit['scale'][0],1.)
        self.assertEqual(fit['mean'][11:],[1.]*3)

    def test_rank_invariance_and_no_support_knots_in_fit_records(self):
        core,proteins,y,w,clusters,folds=synthetic(n=200)
        fit=m.fit_preprocessing(core,proteins,w,'rank_all3')
        other=m.fit_preprocessing(core,np.exp(.2*proteins),w,'rank_all3')
        np.testing.assert_array_equal(m.transform_features(core,proteins,fit),
                                     m.transform_features(core,np.exp(.2*proteins),other))
        _,records=m.nested_crossfit(core,proteins,y,w,clusters,folds,'rank_all3')
        for row in records:
            for diagnostic in row['protein_transform_parameters']:
                self.assertNotIn('values',diagnostic);self.assertNotIn('cdf',diagnostic)
                self.assertGreater(diagnostic['training_support_count'],0)

    def test_legacy_registries_and_inputs_are_not_mutated(self):
        registry=copy.deepcopy(old.MODEL_INDICES);names=old.CORE_NAMES
        core,protein,y,w,clusters,folds=synthetic(n=200)
        originals=[a.copy() for a in (core,protein,y,w,clusters,folds)]
        for name in m.MODEL_INDICES:m.nested_crossfit(core,protein,y,w,clusters,folds,name)
        self.assertEqual(old.MODEL_INDICES,registry);self.assertEqual(old.CORE_NAMES,names)
        for current,original in zip((core,protein,y,w,clusters,folds),originals):np.testing.assert_array_equal(current,original)


class NestedTests(unittest.TestCase):
    def test_simple_predictions_ignore_unselected_routine_features(self):
        core,protein,y,w,clusters,folds=synthetic()
        other=core.copy();other[:,6:]*=1e4
        for name in ('simple','simple_all3'):
            p,a=m.nested_crossfit(core,protein,y,w,clusters,folds,name)
            q,b=m.nested_crossfit(other,protein,y,w,clusters,folds,name)
            np.testing.assert_array_equal(p,q);self.assertEqual(a,b)

    def test_heldout_outcomes_cannot_change_own_prediction_or_selection(self):
        core,protein,y,w,clusters,folds=synthetic()
        p,a=m.nested_crossfit(core,protein,y,w,clusters,folds,'routine_all3')
        changed=y.copy();changed[folds==0]=1-changed[folds==0]
        q,b=m.nested_crossfit(core,protein,changed,w,clusters,folds,'routine_all3')
        np.testing.assert_array_equal(p[folds==0],q[folds==0])
        for key in ('selected_penalty','training_mean','training_scale','inner_brier_by_penalty'):
            self.assertEqual(a[0][key],b[0][key])

    def test_heldout_routine_values_never_enter_training_preprocessing(self):
        core,protein,y,w,clusters,folds=synthetic()
        core[:,0]=np.arange(len(core));calls=[];actual=m.fit_preprocessing
        def capture(c,p,weight,name):
            calls.append(set(c[:,0].astype(int)))
            return actual(c,p,weight,name)
        with patch.object(m,'fit_preprocessing',side_effect=capture):
            m.nested_crossfit(core,protein,y,w,clusters,folds,'routine_all3')
        heldout=set(np.flatnonzero(folds==0))
        self.assertEqual(len(calls),20)
        for call in calls[:4]:self.assertFalse(call&heldout)

    def test_psu_may_not_cross_outer_boundary(self):
        args=list(synthetic());args[-1][1]=(args[-1][1]+1)%5
        with self.assertRaisesRegex(ValueError,'PSU crosses'):m.nested_crossfit(*args)

    def test_all_arms_share_fold_support_and_valid_predictions(self):
        args=synthetic(n=200);supports=[]
        for name in m.MODEL_INDICES:
            p,rows=m.nested_crossfit(*args,model_name=name)
            self.assertTrue(np.isfinite(p).all());self.assertTrue(((p>=0)&(p<=1)).all())
            supports.append([(r['outer_fold'],r['training_n'],r['scoring_n'],r['training_events'],r['scoring_events']) for r in rows])
            expected=6 if name.startswith('simple') else 11
            for row in rows:
                self.assertEqual(row['core_feature_count'],expected)
                self.assertEqual(len(row['core_coefficients_standardized']),expected)
                self.assertEqual(row['PSU_overlap_n'],0);self.assertLessEqual(row['gradient_max'],1e-4)
        self.assertTrue(all(support==supports[0] for support in supports))

    def test_zero_weight_psus_keep_original_fold_assignment(self):
        core,protein,y,w,clusters,folds=synthetic()
        w[np.isin(clusters,sorted(set(clusters))[:10])]=0
        _,rows=m.nested_crossfit(core,protein,y,w,clusters,folds,'routine_all3')
        for row in rows:
            self.assertLess(row['training_positive_weight_n'],row['training_n'])
            self.assertAlmostEqual(row['normalized_weight_sum'],row['training_n'])

    def test_ties_choose_stronger_predeclared_ridge(self):
        core,protein,y,w,clusters,folds=synthetic();core[:]=0;protein[:]=0
        _,rows=m.nested_crossfit(core,protein,y,w,clusters,folds,'routine_all3')
        self.assertTrue(all(r['lambda']==100. for r in rows))

    def test_simple_numerics_match_binary28_with_same_inner_partition(self):
        core,protein,y,w,clusters,folds=synthetic()
        # Force only the inner partition salt to the original for this numerical
        # equivalence check; routine29's production salt is distinct by design.
        original=m.psu_partition
        def old_salt(c,salt,folds=3):return original(c,salt.replace('routine29:','binary28:'),folds)
        for name,legacy in [('simple','core'),('simple_all3','all3')]:
            with patch.object(m,'psu_partition',side_effect=old_salt):
                p,_=m.nested_crossfit(core,protein,y,w,clusters,folds,name)
            q,_=old.nested_crossfit(core[:,:6],protein,y,w,clusters,folds,legacy)
            np.testing.assert_allclose(p,q,atol=1e-14,rtol=0)

if __name__=='__main__':unittest.main()
