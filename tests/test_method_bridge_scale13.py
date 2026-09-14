"""Extra numerical regression for finite, artificially large paired values."""
import unittest
from research.mortality.method_bridge13 import agreement

class LargeScaleTests(unittest.TestCase):
    def test_correlation_large_finite_scale(self):
        result=agreement([1e80,2e80,3e80],[3e80,5e80,7e80])
        self.assertAlmostEqual(result['pearson_correlation'],1.)

if __name__=='__main__':unittest.main()
