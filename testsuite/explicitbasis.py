import unittest

import numpy as np

from pairinteraction import picomplex as pi


class ExplicitBasisTest(unittest.TestCase):
    def setUp(self):
        # Set up cache
        self.cache = pi.MatrixElementCache()

    def test_basisvectors(self):

        one_atom_basisvectors_indices = [[0, 0], [0, 1], [1, 0], [1, 1]]

        # Setup states
        state_one = pi.StateOne("Cs", 60, 0, 0.5, 0.5)

        # One-atom system
        system_one = pi.SystemOne(state_one.getSpecies(), self.cache)
        system_one.restrictEnergy(state_one.getEnergy() - 30, state_one.getEnergy() + 30)
        system_one.restrictN(state_one.getN() - 1, state_one.getN() + 1)
        system_one.restrictL(state_one.getL() - 1, state_one.getL() + 1)
        system_one.setEfield([5, 0, 0])
        system_one.diagonalize()
        basisvectors_one = system_one.getBasisvectors().toarray()
        assert basisvectors_one.shape[0] == basisvectors_one.shape[1]

        # Two-atom system
        system_two = pi.SystemTwo(system_one, system_one, self.cache)
        system_two.setMinimalNorm(0)
        system_two.setOneAtomBasisvectors(one_atom_basisvectors_indices)
        basisvectors_two = system_two.getBasisvectors().toarray()
        assert basisvectors_two.shape[0] == (basisvectors_one.shape[0]) ** 2
        assert basisvectors_two.shape[1] == len(one_atom_basisvectors_indices)

        # Check results
        for n, [a, b] in enumerate(one_atom_basisvectors_indices):
            np.testing.assert_allclose(
                np.outer(basisvectors_one[:, a], basisvectors_one[:, b]).ravel(), basisvectors_two[:, n]
            )

    def test_exception(self):

        one_atom_basisvectors_indices = [[0, 0], [0, 0]]

        # Setup states
        state_one = pi.StateOne("Cs", 60, 0, 0.5, 0.5)

        # One-atom system
        system_one = pi.SystemOne(state_one.getSpecies(), self.cache)
        system_one.restrictEnergy(state_one.getEnergy() - 30, state_one.getEnergy() + 30)
        system_one.restrictN(state_one.getN() - 1, state_one.getN() + 1)
        system_one.restrictL(state_one.getL() - 1, state_one.getL() + 1)
        system_one.setEfield([5, 0, 0])
        system_one.diagonalize()
        _ = system_one.getBasisvectors().toarray()

        # Two-atom system
        system_two = pi.SystemTwo(system_one, system_one, self.cache)
        system_two.setMinimalNorm(0)
        with self.assertRaises(RuntimeError) as context:
            system_two.setOneAtomBasisvectors(one_atom_basisvectors_indices)
        self.assertTrue("not unique" in str(context.exception))


if __name__ == "__main__":
    unittest.main()
