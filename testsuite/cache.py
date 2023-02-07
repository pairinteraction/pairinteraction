import numpy as np
import unittest

from @LIBNAME@ import pireal as pi


class CacheTest(unittest.TestCase):

    def test_defectdb(self):
        cache = pi.MatrixElementCache()
        cache.setDefectDB("defects.sql")

    def test_dipoledb(self):
        cache = pi.MatrixElementCache()
        cache.loadElectricDipoleDB("dipole.csv", "Rb")
        cache_size = cache.size()
        self.assertGreater(cache_size, 0)

        cache_comparison = pi.MatrixElementCache()
        cache_comparison_size = cache_comparison.size()
        self.assertEqual(cache_comparison_size, 0)

        state_f = pi.StateOne("Rb", 8, 0, 1 / 2, 1 / 2)
        state_i = pi.StateOne("Rb", 8, 1, 1 / 2, 1 / 2)
        radial = cache.getRadial(state_f, state_i, 1)
        radial_comparison = cache_comparison.getRadial(state_f, state_i, 1)
        self.assertAlmostEqual(radial, -0.0016741253849814623, places=12)
        self.assertAlmostEqual(radial, radial_comparison, places=5)
        self.assertEqual(cache.size(), cache_size)
        self.assertEqual(cache_comparison.size(), cache_comparison_size + 1)

    @unittest.skipIf(not pi.gsl_enabled, "The program was compiled without GSL support.")
    def test_radial_methods(self):
        cache_numerov = pi.MatrixElementCache()
        cache_numerov.setMethod(pi.NUMEROV)
        cache_whittaker = pi.MatrixElementCache()
        cache_whittaker.setMethod(pi.WHITTAKER)

        state = pi.StateOne("Rb", 42, 0, 1 / 2, 1 / 2)
        for kappa in range(10):
            radial_numerov = cache_numerov.getRadial(state, state, kappa)
            radial_whittaker = cache_whittaker.getRadial(state, state, kappa)
            self.assertAlmostEqual(radial_numerov, radial_whittaker, places=5)

        self.assertEqual(cache_numerov.size(), 10)
        self.assertEqual(cache_whittaker.size(), 10)


if __name__ == '__main__':
    unittest.main()
