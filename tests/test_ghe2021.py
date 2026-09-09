import json, math, os, tempfile, unittest
from pathlib import Path
from unittest.mock import patch
from research.mortality.ghe2021 import (
    assert_antichain, build, check_close, disease_cut, rank_bounds, read_sheet)

BASE=Path(__file__).resolve().parents[1]/'research/mortality'

class RankTests(unittest.TestCase):
    def test_hand_bounds(self):
        r=rank_bounds({1:100,2:80,3:20,4:90,5:30},[1,2,3],[4,5])
        self.assertEqual([(x['rank_lower'],x['rank_upper'])for x in r],[(1,1),(2,3),(3,5)])
    def test_no_mixed(self):
        r=rank_bounds({1:2,2:3},[1,2],[])
        self.assertEqual([x['ghe_code']for x in r],[2,1])
        self.assertTrue(all(x['rank_lower']==x['rank_upper']for x in r))
    def test_tie_rule(self):
        r=rank_bounds({1:10,2:10,3:10},[2,3],[1])
        self.assertEqual([(x['rank_lower'],x['rank_upper'])for x in r],[(1,2),(2,3)])
    def test_invalid_bounds(self):
        for counts,compatible,mixed in [({1:-1},[1],[]),({1:math.nan},[1],[]),
                ({1:1},[2],[]),({1:1},[1],[1]),({1:1},[1,1],[])]:
            with self.subTest(counts=counts,compatible=compatible,mixed=mixed),self.assertRaises(ValueError):
                rank_bounds(counts,compatible,mixed)

class HierarchyTests(unittest.TestCase):
    def setUp(self):
        self.nodes={0:{'depth':-1,'children':[10]},10:{'depth':0,'children':[20,30]},
                    20:{'depth':1,'children':[21,22]},21:{'depth':2,'children':[211]},
                    211:{'depth':3,'children':[]},22:{'depth':2,'children':[]},
                    30:{'depth':1,'children':[]}}
    def test_cut(self):self.assertEqual(set(disease_cut(self.nodes)),{21,22,30})
    def test_double_count_rejected(self):
        with self.assertRaises(ValueError):assert_antichain(self.nodes,[20,21])
    def test_duplicates_rejected(self):
        with self.assertRaises(ValueError):assert_antichain(self.nodes,[21,21])
    def test_unknown_rejected(self):
        with self.assertRaises(ValueError):assert_antichain(self.nodes,[99])
    def test_failed_reconciliation(self):
        with self.assertRaises(ValueError):check_close(100,101,'test')

class ReaderTests(unittest.TestCase):
    def fixture(self):
        return [(8,{'X':'RUS','Y':'DEU','Z':'USA'}),
          (10,{'A':'Persons','D':"Population ('000)",'X':'100','Y':'200','Z':'300'}),
          (11,{'A':'Persons','B':'0','D':'All Causes','X':'0.2','Y':'0.4','Z':'0.6'}),
          (12,{'A':'Persons','B':'1510','C':'III.','D':'Injuries','X':'0.1','Y':'0.2','Z':'0.3'}),
          (13,{'A':'Persons','B':'1700','D':'C.','E':'Other pandemic','X':'0.1','Y':'0.2','Z':'0.3'})]
    def test_thousands_and_dynamic_columns(self):
        with patch('research.mortality.ghe2021.iter_rows',return_value=iter(self.fixture())):
            records,_=read_sheet(Path('dummy'),2,metric='deaths',year=2021,age='All ages')
        r=next(x for x in records if x['country']=='RUS'and x['ghe_code']==0)
        self.assertEqual(r['value'],200);self.assertEqual(r['population'],100000)
        self.assertEqual(r['source_cell'],'X11')
    def test_pandemic_is_not_injuries_child(self):
        with patch('research.mortality.ghe2021.iter_rows',return_value=iter(self.fixture())):
            _,nodes=read_sheet(Path('dummy'),2,metric='deaths',year=2021,age='All ages')
        self.assertEqual(nodes[1700]['parent'],0);self.assertEqual(nodes[1510]['children'],[])
    def test_rate_not_multiplied_by_thousand(self):
        with patch('research.mortality.ghe2021.iter_rows',return_value=iter(self.fixture())):
            records,_=read_sheet(Path('dummy'),2,metric='asdr',year=2021,age='All ages')
        self.assertEqual(records[0]['value'],.2)
    def test_missing_country_rejected(self):
        with patch('research.mortality.ghe2021.iter_rows',return_value=iter([(8,{'X':'RUS'})])):
            with self.assertRaises(ValueError):read_sheet(Path('dummy'),2,metric='deaths',year=2021,age='All ages')

class PolicyTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.policy=json.loads((BASE/'ghe2021_policy.json').read_text())
    def test_partition_policy_size(self):self.assertEqual(len(self.policy['causes']),134)
    def test_external_poisonings_not_compatible(self):
        for c in ('860','870'):
            self.assertEqual(self.policy['causes'][c]['eligibility'],'mixed_or_unreviewed')
            self.assertIn('X4',self.policy['causes'][c]['icd10'])
    def test_no_false_kidney_exclusion(self):
        self.assertEqual(self.policy['causes']['1270']['eligibility'],'mixed_or_unreviewed')
        self.assertIn('E10.2',self.policy['causes']['1270']['icd10'])
    def test_noninfectious_not_identical_to_ncd(self):
        self.assertEqual(self.policy['causes']['550']['eligibility'],'compatible')
        self.assertEqual(self.policy['causes']['520']['eligibility'],'excluded_infection')
    def test_pandemic_unallocated(self):
        self.assertEqual(self.policy['causes']['1700']['eligibility'],'unallocated')

class SourceIntegrationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        path=os.environ.get('GHE_SOURCE_DIR')
        if not path:raise unittest.SkipTest('Set GHE_SOURCE_DIR to downloaded, hash-pinned WHO inputs')
        cls.tmp=tempfile.TemporaryDirectory()
        cls.audit=build(Path(path),BASE/'ghe2021_policy.json',Path(cls.tmp.name))
        cls.data=json.loads((Path(cls.tmp.name)/'rankings.json').read_text())
    @classmethod
    def tearDownClass(cls):cls.tmp.cleanup()
    def test_full_inventory(self):
        self.assertEqual(self.audit['observations'],23220)
        self.assertEqual(self.audit['partition_groups'],134)
    def test_numerical_audits(self):
        self.assertEqual(sum(self.audit['checks'].values()),13026)
        self.assertLess(max(self.audit['max_absolute_errors'].values()),1e-6)
    def test_three_top20(self):
        for iso in ('RUS','DEU','USA'):
            rs=[r for r in self.data['top20_compatible']if r['country']==iso]
            self.assertEqual(len(rs),20);self.assertEqual(len({r['ghe_code']for r in rs}),20)
            self.assertEqual(rs[0]['ghe_code'],1130)
            self.assertTrue(all(r['rank_lower']<=r['rank_upper']for r in rs))
    def test_russian_pandemic_residual_preserved(self):
        r=next(x for x in self.data['countries']if x['country']=='RUS')
        self.assertAlmostEqual(r['unallocated'],161872.76306467,places=5)
    def test_no_individual_probabilities(self):
        self.assertIn('not_individual_risk',self.data['status'])

if __name__=='__main__':unittest.main()
