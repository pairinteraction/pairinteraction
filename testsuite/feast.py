import unittest

import numpy as np

from pairinteraction import picomplex as pi


class FeastTest(unittest.TestCase):
    def setUp(self):
        # Set up cache
        self.cache = pi.MatrixElementCache()

    @unittest.skipIf(not pi.mkl_enabled, "The program was compiled without MKL support.")
    def test_diagonalization_full(self):
        # Setup states
        state_one = pi.StateOne("Cs", 60, 0, 0.5, 0.5)

        # Build one-atom system
        system_one = pi.SystemOne(state_one.getSpecies(), self.cache)
        system_one.restrictEnergy(state_one.getEnergy() - 30, state_one.getEnergy() + 30)
        system_one.restrictN(state_one.getN() - 1, state_one.getN() + 1)
        system_one.restrictL(state_one.getL() - 1, state_one.getL() + 1)
        system_one.setEfield([1, 0, 0])
        system_one.setBfield([5, 50, 0])

        # Diagonalize using the standard approach
        system_one_standard = pi.SystemOne(system_one)
        system_one_standard.diagonalize()
        evecs_standard = system_one_standard.getBasisvectors()

        # Diagonalize using FEAST
        system_one.diagonalize(-1e12, 1e12)
        evecs_feast = system_one.getBasisvectors()

        # Check results
        overlap = np.abs(np.dot(evecs_standard.conj().T, evecs_feast)) ** 2
        np.testing.assert_allclose(overlap.diagonal(), np.ones_like(overlap.diagonal()), rtol=1e-12)
        self.assertAlmostEqual(np.sum(overlap), system_one.getNumBasisvectors(), places=6)

    @unittest.skipIf(not pi.mkl_enabled, "The program was compiled without MKL support.")
    def test_diagonalization_bounded(self):
        # Setup states
        state_one = pi.StateOne("Cs", 60, 0, 0.5, 0.5)

        # Build one-atom system
        system_one = pi.SystemOne(state_one.getSpecies(), self.cache)
        system_one.restrictEnergy(state_one.getEnergy() - 100, state_one.getEnergy() + 100)
        system_one.restrictN(state_one.getN() - 1, state_one.getN() + 1)
        system_one.restrictL(state_one.getL() - 1, state_one.getL() + 1)
        system_one.restrictM(state_one.getM(), state_one.getM())
        system_one.setEfield([0, 0, 1])
        system_one.setBfield([0, 0, 50])

        # Determine energy bounds
        energy_lower_bound = state_one.getEnergy() - 20
        energy_upper_bound = state_one.getEnergy() + 20

        # Diagonalize using FEAST
        system_one.diagonalize(energy_lower_bound, energy_upper_bound, 0.01)

        # Check results
        self.assertEqual(system_one.getNumStates(), 9)
        self.assertEqual(system_one.getNumBasisvectors(), 4)


if __name__ == "__main__":
    unittest.main()
