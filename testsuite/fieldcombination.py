import numpy as np
import unittest

from pairinteraction import pireal as pi


class FieldCombinationTest(unittest.TestCase):

    def test_combined_fields(self):
        cache = pi.MatrixElementCache()

        # Set up SystemOne
        system_one = pi.SystemOne("Rb", cache)
        system_one.restrictEnergy(-1077.243011609127, -939.9554235203701)
        system_one.restrictN(57, 63)
        system_one.restrictL(0, 3)
        system_one.setConservedMomentaUnderRotation([-0.5])
        system_one.setEfield((0, 0, 0.7))
        system_one.setBfield((0, 0, -8.8))
        system_one.enableDiamagnetism(False)

        # Diagonalize the system
        system_one.diagonalize()

        # Compare results
        energies = system_one.getHamiltonian().diagonal()
        self.assertAlmostEqual(energies[13], -1000.2679341660352, places=4)


if __name__ == '__main__':
    unittest.main()
