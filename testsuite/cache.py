import numpy as np
import pickle
import shutil
import tempfile
import unittest
import ssl
import urllib.request
import os

from @LIBNAME@ import pireal as pi


class CacheTest(unittest.TestCase):

    def setUp(self):
        self.cache_path = tempfile.mkdtemp()
        self.dipoledb_path = os.path.join(self.cache_path, "dipole.csv")
        self.defectdb_path = os.path.join(self.cache_path, "defects.sql")

        ctx = ssl.create_default_context()
        ctx.check_hostname = False
        ctx.verify_mode = ssl.CERT_NONE

        with urllib.request.urlopen("https://raw.githubusercontent.com/nikolasibalic/ARC-Alkali-Rydberg-Calculator/4579309f65152fa65ca4fc5796eb00d199a4dc39/arc/data/rubidium_literature_dme.csv", context=ctx) as u, \
                open(self.dipoledb_path, 'wb') as f:
            f.write(u.read())
        with urllib.request.urlopen("https://raw.githubusercontent.com/pairinteraction/pairinteraction/dcd4013238984ea2f57ffc2cd1eaff2996fc8586/libpairinteraction/databases/quantum_defects.sql", context=ctx) as u, \
                open(self.defectdb_path, 'wb') as f:
            f.write(u.read())

    def tearDown(self):
        try:
            shutil.rmtree(self.cache_path)
        except BaseException:
            pass

    def test_defectdb(self):
        cache = pi.MatrixElementCache()
        cache.setDefectDB(self.defectdb_path)

    def test_dipoledb(self):
        cache = pi.MatrixElementCache()
        cache.loadElectricDipoleDB(self.dipoledb_path, "Rb")
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
