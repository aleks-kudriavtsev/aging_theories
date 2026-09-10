"""Independent arithmetic and adversarial controls; artificial pairs only."""
from copy import deepcopy
import math
import tempfile
from pathlib import Path
import unittest
from research.mortality.method_bridge13 import (deming,correct,agreement,analyse,demo,
    validate_document,application_gate,numeric,quantile,load_pairs_csv)

class DemingTests(unittest.TestCase):
    def test_affine_exact(self):
        m=deming([1,2,3,4],[5,7,9,11],variance_ratio=3)
        self.assertAlmostEqual(m['slope'],2);self.assertAlmostEqual(m['intercept'],3)
    def test_variance_ratio_is_required(self):
        with self.assertRaises(TypeError):deming([1,2,3],[2,3,4])
    def test_inverse_fit_consistency(self):
        x=[1,3,4,7,9];y=[2,4,5,10,14]
        a=deming(x,y,variance_ratio=4);b=deming(y,x,variance_ratio=.25)
        self.assertAlmostEqual(a['slope']*b['slope'],1)
        self.assertAlmostEqual(b['intercept'],-a['intercept']/a['slope'])
    def test_reject_constant(self):
        with self.assertRaises(ValueError):deming([1,1,1],[2,3,4],variance_ratio=1)
    def test_reject_negative_relationship(self):
        with self.assertRaises(ValueError):deming([1,2,3],[3,2,1],variance_ratio=1)
    def test_reject_mismatched_arrays(self):
        with self.assertRaises(ValueError):deming([1,2,3],[3,2],variance_ratio=1)
    def test_reject_nonnumeric(self):
        for v in [True,math.nan,math.inf,None,'<3',10**1000,0,-1]:
            with self.subTest(v=str(v)[:20]),self.assertRaises(ValueError):deming([1,2,v],[2,3,4],variance_ratio=1)
    def test_ratio_bounds(self):
        for v in [None,True,0,-1,math.nan,1e100]:
            with self.subTest(v=v),self.assertRaises(ValueError):deming([1,2,3],[2,3,4],variance_ratio=v)
    def test_correction_exact(self):
        m=deming([1,2,3],[3,5,7],variance_ratio=1)
        self.assertAlmostEqual(correct(5,m),2)
    def test_correction_extrapolation_forbidden(self):
        m=deming([1,2,3],[3,5,7],variance_ratio=1)
        with self.assertRaises(ValueError):correct(7.1,m)
    def test_invalid_correction_spec(self):
        for m in [None,{}, {'slope':0,'intercept':1,'candidate_range':[1,2]}]:
            with self.subTest(m=m),self.assertRaises(ValueError):correct(1,m)
    def test_limits_of_agreement_are_not_ci(self):
        r=agreement([1,2,3],[2,3,4])
        self.assertEqual(r['mean_difference'],1);self.assertEqual(r['approximate_limits_of_agreement'],[1,1])
        self.assertIn('NOT confidence',r['limits_definition'])
    def test_perfect_correlation_can_have_bias(self):
        r=agreement([1,2,3],[3,5,7])
        self.assertAlmostEqual(r['pearson_correlation'],1);self.assertEqual(r['mean_difference'],3)
        self.assertFalse(r['correlation_establishes_interchangeability'])
    def test_numeric_overflow_handled(self):self.assertFalse(numeric(10**1000))
    def test_quantile(self):self.assertEqual(quantile([1,3],.5),2)

class DocumentTests(unittest.TestCase):
    def block(self,mutate):
        d=demo();mutate(d)
        with self.assertRaises(ValueError):validate_document(d)
    def test_demo_splits(self):
        s=validate_document(demo());self.assertEqual(len(s['fit']),61);self.assertEqual(len(s['validation']),46)
    def test_duplicate_specimen(self):self.block(lambda d:d['pairs'].append(deepcopy(d['pairs'][0])))
    def test_cross_split_leak(self):
        def mutate(d):d['pairs'][-1]['specimen_id']=d['pairs'][0]['specimen_id']
        self.block(mutate)
    def test_invalid_split_mapping(self):
        def mutate(d):d['pairs'][0]['split']=[]
        self.block(mutate)
    def test_specimen_name_extra_rejected(self):self.block(lambda d:d['pairs'][0].update(patient_name='example'))
    def test_extra_top_level(self):self.block(lambda d:d.update(outcome='death'))
    def test_unknown_origin(self):self.block(lambda d:d.update(data_origin='real_verified'))
    def test_different_analyte(self):self.block(lambda d:d['candidate_method'].update(analyte='other'))
    def test_different_matrix(self):self.block(lambda d:d['candidate_method'].update(matrix='other'))
    def test_different_unit(self):self.block(lambda d:d['candidate_method'].update(unit='other'))
    def test_unstated_lot(self):self.block(lambda d:d['candidate_method'].update(lot_id='not_reported'))
    def test_variance_ratio_basis(self):self.block(lambda d:d.update(variance_ratio_basis=''))
    def test_unknown_risk_model(self):self.block(lambda d:d.update(model_sha256=''))
    def test_sample_size_not_automatic_validity(self):self.block(lambda d:d.update(pairs=d['pairs'][:5]))
    def test_unquantified_result_not_zero(self):
        def mutate(d):d['pairs'][0]['candidate']='<2'
        self.block(mutate)

class AnalysisTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.doc=demo();cls.result=analyse(cls.doc,bootstrap_replicates=100)
    def test_synthetic_cannot_become_validation(self):
        r=self.result;self.assertEqual(r['status'],'synthetic_method_demo');self.assertFalse(r['mortality_transport_validated'])
        self.assertFalse(r['clinical_use_ready']);self.assertFalse(r['automatic_risk_model_activation'])
    def test_known_synthetic_improvement(self):
        r=self.result;self.assertLess(r['validation_corrected']['mean_absolute_difference'],r['validation_raw']['mean_absolute_difference'])
    def test_no_specimen_rows_exported(self):
        self.assertNotIn('pairs',self.result);self.assertFalse(self.result['specimen_rows_exported'])
    def test_validation_does_not_change_fit(self):
        d=demo()
        for row in d['pairs']:
            if row['split']=='validation':row['candidate']+=.01
        r=analyse(d,bootstrap_replicates=100);self.assertEqual(r['fit'],self.result['fit'])
        self.assertNotEqual(r['validation_raw'],self.result['validation_raw'])
    def test_outside_range_not_silently_excluded(self):
        d=demo();d['pairs'][-1]['candidate']=100
        r=analyse(d,bootstrap_replicates=100)
        self.assertEqual(r['validation_outside_fitted_range'],1);self.assertIsNone(r['validation_corrected'])
        self.assertEqual(r['validation_raw']['n'],46)
    def test_no_mutation(self):
        d=demo();saved=deepcopy(d);analyse(d,bootstrap_replicates=100);self.assertEqual(d,saved)
    def test_bootstrap_reproducible(self):self.assertEqual(analyse(demo(),bootstrap_replicates=100),self.result)
    def test_bootstrap_bounds(self):
        for v in [True,0,99,10001]:
            with self.subTest(v=v),self.assertRaises(ValueError):analyse(demo(),bootstrap_replicates=v)
    def test_no_automatic_application_even_if_forged_flag(self):
        r=deepcopy(self.result);r['clinical_use_ready']=True
        m=self.doc['candidate_method'];g=application_gate(r,**m,model_sha256=self.doc['model_sha256'])
        self.assertFalse(g['ready']);self.assertIn('synthetic_or_undeclared_data',g['reasons'])
    def test_application_context_checked(self):
        m=deepcopy(self.doc['candidate_method']);m['matrix']='different'
        g=application_gate(self.result,**m,model_sha256='0'*64)
        self.assertIn('different_matrix',g['reasons']);self.assertIn('different_risk_model',g['reasons'])
    def test_actual_declaration_is_not_source_verification(self):
        d=demo();d['data_origin']='paired_specimens_declared';r=analyse(d,bootstrap_replicates=100)
        self.assertFalse(r['source_independently_verified']);self.assertFalse(r['mortality_transport_validated'])

class CsvTests(unittest.TestCase):
    def check_csv(self,text):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'pairs.csv';p.write_text(text,encoding='utf-8');return load_pairs_csv(p)
    def test_csv_read(self):
        r=self.check_csv('specimen_id,split,reference,candidate\nS1,fit,1,2\n')
        self.assertEqual(r[0]['reference'],1.)
    def test_csv_wrong_header(self):
        with self.assertRaises(ValueError):self.check_csv('patient_id,split,reference,candidate\n')
    def test_csv_unquantified(self):
        for v in ['', 'NaN', 'inf', '<3']:
            with self.subTest(v=v),self.assertRaises(ValueError):self.check_csv('specimen_id,split,reference,candidate\nS1,fit,1,'+v+'\n')
    def test_csv_extra_column(self):
        with self.assertRaises(ValueError):self.check_csv('specimen_id,split,reference,candidate\nS1,fit,1,2,extra\n')
    def test_csv_missing_column(self):
        with self.assertRaises(ValueError):self.check_csv('specimen_id,split,reference,candidate\nS1,fit,1\n')

if __name__=='__main__':unittest.main()
