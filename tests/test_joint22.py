"""Artificial fixtures only; no HRS or UKB participants."""
from copy import deepcopy
import unittest
from research.mortality.nfl22.joint_data import PROTEINS,validate_plan,audit_overlap

def fixture():
    methods={m:{'source_variable':'fixture_'+m,'source_analyte':m,'identity_verified':True,
                'matrix':'plasma','unit':'NPX','scale':'NPX','assay':'synthetic_assay'} for m in PROTEINS}
    plan={'cohort':'SYNTHETIC','baseline_visit':'v0','approval_reference':'test_only_not_permission',
          'data_release':'test','dictionary_sha256':'test_hash','data_use_approved':True,
          'required_markers':list(PROTEINS),'measurements':methods}
    rows=[{'cohort':'SYNTHETIC','participant_id':'test_1','visit':'v0','marker':m,
           'matrix':'plasma','unit':'NPX','assay':'synthetic_assay','qualifier':'quantified','value':.5} for m in PROTEINS]
    return rows,plan

class JoinTests(unittest.TestCase):
    def test_exact_overlap(self):
        rows,p=fixture();r=audit_overlap(rows,p);self.assertEqual(r['complete_quantified_profiles'],1);self.assertFalse(r['identifiers_exported'])
    def test_no_approval(self):
        _,p=fixture();p['data_use_approved']=False
        with self.assertRaises(ValueError):validate_plan(p)
    def test_receptor_is_not_ligand(self):
        _,p=fixture();p['measurements']['TNF_ALPHA']['source_analyte']='TNFRSF1A'
        with self.assertRaises(ValueError):validate_plan(p)
    def test_NPPB_not_NT_fragment(self):
        _,p=fixture();p['measurements']['NT_PROBNP']['source_analyte']='NPPB'
        with self.assertRaises(ValueError):validate_plan(p)
    def test_NPX_not_pg_ml(self):
        _,p=fixture();p['measurements']['GFAP']['unit']='pg/mL'
        with self.assertRaises(ValueError):validate_plan(p)
    def test_no_cross_participant_fill(self):
        rows,p=fixture();rows[0]['participant_id']='test_2'
        self.assertEqual(audit_overlap(rows,p)['complete_quantified_profiles'],0)
    def test_no_cross_cohort(self):
        rows,p=fixture();rows[0]['cohort']='OTHER'
        with self.assertRaises(ValueError):audit_overlap(rows,p)
    def test_no_cross_visit(self):
        rows,p=fixture();rows[0]['visit']='v1'
        with self.assertRaises(ValueError):audit_overlap(rows,p)
    def test_censored_retained_not_normal(self):
        rows,p=fixture();rows[0].update(qualifier='below_LOQ',value=None)
        r=audit_overlap(rows,p);self.assertEqual(r['complete_quantified_profiles'],0);self.assertEqual(r['qualifier_counts']['below_LOQ'],1)
    def test_duplicate_not_silent_average(self):
        rows,p=fixture();rows.append(deepcopy(rows[0]))
        with self.assertRaises(ValueError):audit_overlap(rows,p)
    def test_method_mismatch(self):
        rows,p=fixture();rows[0]['matrix']='serum'
        with self.assertRaises(ValueError):audit_overlap(rows,p)
    def test_incomplete_identity_blocked(self):
        rows,p=fixture();p['measurements']['GFAP']['identity_verified']=False
        with self.assertRaises(ValueError):audit_overlap(rows,p)

if __name__=='__main__':unittest.main()
