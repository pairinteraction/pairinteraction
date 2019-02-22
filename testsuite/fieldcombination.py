import numpy as np
import pickle
import shutil
import tempfile
import unittest

from @LIBNAME@ import picomplex as pi


class FieldCombinationTest(unittest.TestCase):

    def test_combined_fields(self):
        cache = pi.MatrixElementCache()

        # Set up SystemOne
        system_one = pi.SystemOne("Rb", cache)
        system_one.restrictEnergy(-1077.243011609127, -939.9554235203701)
        system_one.restrictN(57, 63)
        system_one.restrictL(-2, 3)
        system_one.setEfield((0, 0, 0.7))
        system_one.setBfield((0, 0, -8.8))

        # Diagonalize the system
        system_one.diagonalize()

        # Compare results
        energies = np.real(system_one.getHamiltonian().diagonal()[[34, 54, 53]])
        np.testing.assert_array_almost_equal(energies, [-1017.29997025, -1000.2551621, -1000.26792709], decimal=4)


if __name__ == '__main__':
    unittest.main()
