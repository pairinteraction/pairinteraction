import pickle
import shutil
import tempfile
import unittest

import numpy as np

from pairinteraction import picomplex as pi


class IntegrationTest(unittest.TestCase):
    def setUp(self):
        self.dump_new_reference_data = False
        self.tolerance = 1e-6

        # Set up cache
        self.cache_path = tempfile.mkdtemp()
        self.cache = pi.MatrixElementCache(self.cache_path)

        # Setup states
        self.state_one = pi.StateOne("Rb", 61, 2, 1.5, 1.5)
        self.state_two = pi.StateTwo(self.state_one, self.state_one)

    def tearDown(self):
        try:
            shutil.rmtree(self.cache_path)
        except BaseException:
            pass

        if self.dump_new_reference_data:
            with open("integration_test_referencedata.pickle", "wb") as f:
                pickle.dump((self.hamiltonian_one, self.basis_one, self.hamiltonian_two, self.basis_two), f)

        self.assertFalse(self.dump_new_reference_data)

    def test_integration(self):
        # Build one-atom system
        system_one = pi.SystemOne(self.state_one.getSpecies(), self.cache)
        system_one.restrictEnergy(self.state_one.getEnergy() - 40, self.state_one.getEnergy() + 40)
        system_one.restrictN(self.state_one.getN() - 1, self.state_one.getN() + 1)
        system_one.restrictL(self.state_one.getL() - 1, self.state_one.getL() + 1)
        system_one.setEfield([0, 0, 0.1])
        system_one.setBfield([0, 0, 1])
        system_one.enableDiamagnetism(False)

        # Check for correct dimensions
        self.assertEqual(system_one.getNumBasisvectors(), 64)
        self.assertEqual(system_one.getNumStates(), 64)

        # Compare current results to the reference data (the results have to be
        # compared before diagonalization as the order of the eigenvectors is not
        # fixed)
        hamiltonian_one = system_one.getHamiltonian()
        basis_one = system_one.getBasisvectors()
        # without pruning, max_diff_hamiltonian might be infinity due to
        # division by zero
        hamiltonian_one.data *= abs(hamiltonian_one).data > self.tolerance
        hamiltonian_one.eliminate_zeros()
        basis_one.data *= abs(basis_one).data > self.tolerance
        basis_one.eliminate_zeros()

        if self.dump_new_reference_data:
            self.hamiltonian_one = hamiltonian_one.copy()
            self.basis_one = basis_one.copy()
        else:
            with open("integration_test_referencedata.pickle", "rb") as f:
                hamiltonian_one_reference, basis_one_reference, _, _ = pickle.load(f)
                np.testing.assert_allclose(hamiltonian_one.A, hamiltonian_one_reference.A, rtol=self.tolerance)
                np.testing.assert_allclose(basis_one.A, basis_one_reference.A, rtol=self.tolerance)

        # Diagonalize one-atom system
        system_one.diagonalize()

        # Build one-atom system (for this test, system_one has to be diagonal by
        # itself because diagonalization can lead to different order of
        # eigenvectors)
        system_one = pi.SystemOne(self.state_one.getSpecies(), self.cache)
        system_one.restrictEnergy(self.state_one.getEnergy() - 40, self.state_one.getEnergy() + 40)
        system_one.restrictN(self.state_one.getN() - 1, self.state_one.getN() + 1)
        system_one.restrictL(self.state_one.getL() - 1, self.state_one.getL() + 1)

        # Build two-atom system
        system_two = pi.SystemTwo(system_one, system_one, self.cache)
        system_two.restrictEnergy(self.state_two.getEnergy() - 2, self.state_two.getEnergy() + 2)
        system_two.setConservedParityUnderPermutation(pi.ODD)
        system_two.setDistance(6)
        system_two.setAngle(0.9)

        # Check for correct dimensions
        self.assertEqual(system_two.getNumBasisvectors(), 239)
        self.assertEqual(system_two.getNumStates(), 468)

        # Compare current results to the reference data (the results have to be
        # compared before diagonalization as the order of the eigenvectors is not
        # fixed)
        hamiltonian_two = system_two.getHamiltonian()
        basis_two = system_two.getBasisvectors()
        # without pruning, max_diff_hamiltonian might be infinity due to
        # division by zero
        hamiltonian_two.data *= abs(hamiltonian_two).data > self.tolerance
        hamiltonian_two.eliminate_zeros()
        basis_two.data *= abs(basis_two).data > self.tolerance
        basis_two.eliminate_zeros()

        if self.dump_new_reference_data:
            self.hamiltonian_two = hamiltonian_two.copy()
            self.basis_two = basis_two.copy()
        else:
            with open("integration_test_referencedata.pickle", "rb") as f:
                _, _, hamiltonian_two_reference, basis_two_reference = pickle.load(f)
                np.testing.assert_allclose(hamiltonian_two.A, hamiltonian_two_reference.A, rtol=self.tolerance)
                np.testing.assert_allclose(basis_two.A, basis_two_reference.A, rtol=self.tolerance)

        # Diagonalize two-atom system
        system_two.diagonalize()


if __name__ == "__main__":
    unittest.main()
