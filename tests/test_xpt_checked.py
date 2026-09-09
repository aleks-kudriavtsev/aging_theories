"""Byte-level regression for a verified XPORT zero-decoding defect."""
from pathlib import Path
import os
import unittest
import numpy as np
from research.mortality.xpt_checked import exact_zero_rows,read_xpt_checked

class ZeroByteTests(unittest.TestCase):
    def test_positive_zero_exact(self):
        self.assertTrue(exact_zero_rows(bytes(8),0,8,1,0,8)[0])
    def test_signed_zero_exact(self):
        self.assertTrue(exact_zero_rows(bytes.fromhex('8000000000000000'),0,8,1,0,8)[0])
    def test_genuinely_tiny_nonzero_not_clipped(self):
        self.assertFalse(exact_zero_rows(bytes.fromhex('0010000000000000'),0,8,1,0,8)[0])
    def test_missing_dot_not_zero(self):
        self.assertFalse(exact_zero_rows(b'.'+bytes(7),0,8,1,0,8)[0])
    def test_missing_A_not_zero(self):
        self.assertFalse(exact_zero_rows(b'A'+bytes(7),0,8,1,0,8)[0])
    def test_record_stride(self):
        raw=b'header00'+bytes(8)+b'xxxx'+bytes.fromhex('4110000000000000')+b'xxxx'
        np.testing.assert_array_equal(exact_zero_rows(raw,8,12,2,0,8),[True,False])
    def test_truncated_layout_blocked(self):
        with self.assertRaises(ValueError):exact_zero_rows(bytes(8),0,8,2,0,8)
    def test_wrong_numeric_width_blocked(self):
        with self.assertRaises(ValueError):exact_zero_rows(bytes(8),0,8,1,0,1)

class SourceZeroTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        src=os.environ.get('NHANES_SOURCE_DIR')
        if not src:raise unittest.SkipTest('Official NHANES inputs required')
        cls.source=Path(src)
    def test_all26_troponin_zeros_corrected(self):
        d=read_xpt_checked(self.source/'1999_SSTROP_A.xpt')
        self.assertEqual(int(d.SSTNT.eq(0).sum()),26)
        self.assertEqual(int(d.SSTNT.between(0,1e-50,inclusive='neither').sum()),0)
    def test_correction_audit_has_no_identifiers(self):
        d=read_xpt_checked(self.source/'1999_SSTROP_A.xpt');a=d.attrs['xpt_zero_audit']
        self.assertIs(a['epsilon_threshold_used'],False)
        self.assertNotIn('participant_ids',a)
    def test_corrected_cohort_counts(self):
        from research.mortality.nhanes_cohort06 import make_cohort_checked
        d,flow=make_cohort_checked(self.source)
        self.assertEqual(len(d[d.cycle<2003]),4549);self.assertEqual(len(d[d.cycle==2003]),2381)
        self.assertTrue(d.SSTNT.gt(0).all());self.assertGreater(d.SSTNT.min(),1e-50)
        self.assertEqual(int(((d.cycle<2003)&d.dead.eq(1)&d.time_years.le(10)).sum()),756)
        self.assertEqual(int(((d.cycle==2003)&d.dead.eq(1)&d.time_years.le(10)).sum()),395)

if __name__=='__main__':unittest.main()
