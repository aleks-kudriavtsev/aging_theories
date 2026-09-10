"""Synthetic coding controls plus optional real-source arithmetic reconciliation."""
import csv,json,os,tempfile,unittest
from copy import deepcopy
from pathlib import Path
from research.mortality.registered_causes10 import (AGES,POLICY,age_vector,harmonize,map_detailed,
 integer,ranks_and_bounds,endpoint_relation,code_set,validate_selection,build,read_sources)


def row(code='I219',sex='1',n=10,fmt='01',listing='104'):
    d={'Country':'4085','Admin1':'','SubDiv':'','Year':'2019','List':listing,'Cause':code,'Sex':sex,'Frmat':fmt}
    d.update({f'Deaths{i}':'0' for i in range(1,27)});d['Deaths1']=str(n);d['Deaths23']=str(n)
    if fmt=='01':d['Deaths24']=d['Deaths25']=''
    return d


class CodingTests(unittest.TestCase):
    def test_named_causes(self):
        for code,key in [('I219','1067'),('I251','1067'),('I64','1069'),('C349','1034'),('G309','1060'),('J449','1076'),('N189','1086'),('C61','1040')]:
            with self.subTest(code=code):self.assertEqual(map_detailed(code),key)
    def test_infection_and_external(self):
        for code in ['A419','B99']:self.assertEqual(map_detailed(code),'1001')
        for code in ['X45','X429','Y349','V019']:self.assertEqual(map_detailed(code),'1095')
        self.assertEqual(map_detailed('J189'),'1074')
    def test_etiology_not_redefined(self):self.assertEqual(POLICY[map_detailed('C539')]['status'],'compatible')
    def test_malnutrition_not_infection(self):self.assertEqual(POLICY[map_detailed('E46')]['status'],'compatible')
    def test_substance_disorder_not_external_code(self):self.assertEqual(map_detailed('F102'),'1056')
    def test_mixed_group_preserved(self):
        self.assertEqual(map_detailed('I330'),'1068');self.assertEqual(POLICY['1068']['status'],'mixed_or_unresolved')
    def test_nervous_meningitis_remainder(self):self.assertEqual(map_detailed('G019'),'1061')
    def test_condensed_not_detailed(self):
        for value in ['1067','AAA','','i21','I20-25',None]:
            with self.subTest(value=value),self.assertRaises(ValueError):map_detailed(value)
    def test_updated_residual_code(self):self.assertEqual(map_detailed('F03'),'1057');self.assertEqual(map_detailed('I159'),'1071')
    def test_special_codes(self):
        self.assertEqual(map_detailed('U071'),'1903');self.assertEqual(map_detailed('U049'),'1901');self.assertEqual(map_detailed('U099'),'other_special')
    def test_policy_unique_partition(self):self.assertEqual(len(POLICY),65)
    def test_no_invented_COPD_equivalence(self):self.assertIn('proper_subset',endpoint_relation(code_set('J40-J44'),code_set('J40-J47')))
    def test_heart_not_ihd(self):self.assertIn('proper_subset',endpoint_relation(code_set('I20-I25'),code_set('I00-I09','I11','I13','I20-I51')))
    def test_missing_scope_not_equal(self):self.assertEqual(endpoint_relation([],[]),'unresolved')
    def test_equal_scope(self):self.assertEqual(endpoint_relation(['I20'],['I20']),'equal_code_set')


class AgeAndSelectionTests(unittest.TestCase):
    def test_85plus_not85_89(self):
        v=age_vector(row(n=15));self.assertEqual(v[AGES.index('85+')],15)
    def test_merge_older_bands(self):
        r=row(n=15,fmt='00');r.update(Deaths23='8',Deaths24='5',Deaths25='2');self.assertEqual(age_vector(r)[-2],15)
    def test_unknown_age_retained(self):
        r=row(n=15);r.update(Deaths23='10',Deaths26='5');self.assertEqual(age_vector(r)[-1],5)
    def test_infant_detail_not_counted_twice(self):
        r=row(n=5);r.update(Deaths23='0',Deaths2='5',IM_Deaths1='5');self.assertEqual(sum(age_vector(r)),5)
    def test_unsupported_age_format(self):
        with self.assertRaises(ValueError):age_vector(row(fmt='04'))
    def test_bad_age_sum(self):
        r=row();r['Deaths2']='1'
        with self.assertRaises(ValueError):age_vector(r)
    def test_overlapping85plus(self):
        r=row();r['Deaths24']='1'
        with self.assertRaises(ValueError):age_vector(r)
    def test_missing_total_not_zero(self):
        r=row();r['Deaths1']=''
        with self.assertRaises(ValueError):age_vector(r)
    def test_noninteger(self):
        for v in ['1.5','-1','',None]:
            with self.subTest(v=v),self.assertRaises(ValueError):integer(v)
    def test_region_rejected(self):
        r=row();r['SubDiv']='urban'
        with self.assertRaises(ValueError):validate_selection([r])
    def test_year_mixture(self):
        a,b=row(),row();b['Year']='2020'
        with self.assertRaises(ValueError):validate_selection([a,b])
    def test_duplicate(self):
        with self.assertRaises(ValueError):validate_selection([row(),row()])
    def test_overlapping_icd_parent_child(self):
        with self.assertRaises(ValueError):validate_selection([row('I21'),row('I219')])
    def test_unknown_sex_retained(self):
        rows=[row('AAA','1'),row('I219','1'),row('AAA','2'),row('I219','2'),row('AAA','9',1),row('I219','9',1)]
        g,t,a=harmonize(rows);self.assertEqual(sum(t['9']),1)
    def test_partition_sum(self):
        g,t,a=harmonize([row('AAA','1'),row('I219','1'),row('AAA','2'),row('I219','2')]);self.assertEqual(sum(g['1']['1067']),10)
    def test_missing_cause_blocks(self):
        with self.assertRaises(ValueError):harmonize([row('AAA','1'),row('AAA','2')])
    def test_absent_sex_total_blocks(self):
        with self.assertRaises(ValueError):harmonize([row('AAA'),row('I219')])


class BoundsTests(unittest.TestCase):
    def test_uncertainty_not_suppression(self):
        r=ranks_and_bounds({'1067':100,'1069':50,'1068':70});self.assertEqual(r['1069'],{'rank_identified':2,'rank_upper_at_fixed_grouping':3})
    def test_external_cannot_be_target(self):self.assertEqual(ranks_and_bounds({'1067':100,'1095':200})['1067']['rank_upper_at_fixed_grouping'],1)
    def test_ties_competition_rank(self):
        r=ranks_and_bounds({'1067':50,'1069':50});self.assertEqual(r['1069']['rank_identified'],1)
    def test_unknown_can_change_upper_bound(self):self.assertEqual(ranks_and_bounds({'1067':50,'1094':60})['1067']['rank_upper_at_fixed_grouping'],2)


class ActualSourceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        s=os.environ.get('REGISTERED10_SOURCE_DIR')
        if not s:raise unittest.SkipTest('Official registered-mortality sources not provided')
        cls.source=Path(s);cls.temp=tempfile.TemporaryDirectory();cls.out=Path(cls.temp.name)/'output'
        cls.result=build(cls.source,cls.out)
    @classmethod
    def tearDownClass(cls):cls.temp.cleanup()
    def test_common_year_and_rows(self):
        self.assertEqual(self.result['year'],2019);self.assertEqual(self.result['top20_rows'],60);self.assertEqual(self.result['sex_age_rows'],11115)
    def test_national_sum(self):
        rows=list(csv.DictReader((self.out/'country_summary.csv').open()));tot={r['country']:int(r['all_deaths']) for r in rows if r['sex']=='all'}
        self.assertEqual(tot,{'RUS':1798307,'DEU':939520,'USA':2854838})
    def test_status_sums(self):
        for r in csv.DictReader((self.out/'country_summary.csv').open()):
            self.assertEqual(int(r['all_deaths']),sum(int(r[k]) for k in ['compatible','infection','external','ill_defined','mixed_or_unresolved']))
    def test_ru_Alzheimer_identified_not_missing(self):
        r=next(r for r in csv.DictReader((self.out/'detail_identifiability.csv').open()) if r['country']=='RUS' and r['component']=='Alzheimer_G30');self.assertEqual(int(r['deaths']),2302)
    def test_population_discrepancy_not_hidden(self):
        r=[x['age_sum_minus_reported_total'] for x in self.result['population_age_sum_discrepancies'] if x['country']=='DEU'];self.assertEqual(r,[6,6])
    def test_us_age_denominators_not_imputed(self):
        r=[r for r in csv.DictReader((self.out/'sex_age_counts.csv').open()) if r['country']=='USA'];self.assertTrue(all(x['population']=='' for x in r))
    def test_frozen_models_not_fitted(self):self.assertFalse(self.result['new_model_fit']);self.assertFalse(self.result['individual_predictions_changed'])
    def test_no_overwrite(self):
        with self.assertRaises(ValueError):build(self.source,self.out)

if __name__=='__main__':unittest.main()
