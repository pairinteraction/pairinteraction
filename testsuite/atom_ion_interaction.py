import unittest

from pairinteraction import pireal as pi


class AtomIonInteractionTest(unittest.TestCase):
    def test_atom_ion_monopole_dipole(self):
        cache = pi.MatrixElementCache()

        # Set up SystemOneReal
        state = pi.StateOne("Rb", 45, 1, 1.5, 0.5)
        system = pi.SystemOneReal(state.getSpecies(), cache)
        system = pi.SystemOneReal(state.getSpecies(), cache)
        system.restrictEnergy(state.getEnergy() - 100, state.getEnergy() + 100)
        system.restrictN(42, 48)
        system.restrictM(0.5, 0.5)
        system.setIonCharge(1)
        system.setRydIonOrder(1)
        system.setRydIonDistance(1.35)

        # Diagonalize the system
        system.diagonalize()

        # Compare results
        energies = system.getHamiltonian().diagonal() / 29.9792458
        index = system.getBasisvectorIndex(state)
        self.assertAlmostEqual(energies[index], -61.30507832465593, places=4)


if __name__ == "__main__":
    unittest.main()
