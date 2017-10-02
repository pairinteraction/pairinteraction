import unittest
from @LIBNAME@ import pireal as pi
import numpy as np

class TestWavefunction(unittest.TestCase):

    def test_comparison(self):
        qd = pi.QuantumDefect("Rb", 80, 1, 0.5)
        n = pi.Numerov(qd).integrate()
        w = pi.Whittaker(qd).integrate()

        diff = np.sqrt(n[:,0])*n[:,1]-w[:,1]

        for p in diff:
            self.assertAlmostEqual(p, 0, places=3)

if __name__ == '__main__':
    unittest.main()
