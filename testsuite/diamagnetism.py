import numpy as np
import unittest

from @LIBNAME@ import picomplex as pi


class DiamagnetismTest(unittest.TestCase):
    def test_diamagnetism_isotropy(self):
        cache = pi.MatrixElementCache()

        # Set up SystemOne
        state_one = pi.StateOne("Rb", 33, 0, 0.5, 0.5)
        system_one = pi.SystemOne(state_one.getSpecies(), cache)
        system_one.restrictEnergy(
            state_one.getEnergy() - 200, state_one.getEnergy() + 200
        )
        system_one.restrictN(state_one.getN() - 2, state_one.getN() + 2)
        system_one.restrictL(state_one.getL() - 2, state_one.getL() + 2)
        system_one.enableDiamagnetism(True)

        # Calculate energies for manetic fields in different directions
        system_one.setBfield((7e3, 0, 0))
        system_one.diagonalize()
        energies_x = system_one.getHamiltonian().diagonal()

        system_one.setBfield((0, 7e3, 0))
        system_one.diagonalize()
        energies_y = system_one.getHamiltonian().diagonal()

        system_one.setBfield((0, 0, 7e3))
        system_one.diagonalize()
        energies_z = system_one.getHamiltonian().diagonal()

        # Compare results
        np.testing.assert_allclose(energies_x, energies_y, rtol=1e-6)
        np.testing.assert_allclose(energies_x, energies_z, rtol=1e-6)

        # Check the energy shift of state_one
        i1 = system_one.getBasisvectorIndex(pi.StateOne("Rb", 33, 0, 0.5, 0.5))
        i2 = system_one.getBasisvectorIndex(pi.StateOne("Rb", 33, 0, 0.5, -0.5))
        self.assertAlmostEqual(
            energies_z[i1] - state_one.getEnergy(), 19.729554603372435, places=4
        )
        self.assertAlmostEqual(
            energies_z[i2] - state_one.getEnergy(), 0.1120907536969753, places=4
        )


if __name__ == "__main__":
    unittest.main()
