import numpy as np
import pickle
import shutil
import tempfile
import unittest

from @LIBNAME@ import picomplex as pi

class IntegrationTest(unittest.TestCase):

    def setUp(self):
        self.dump_new_reference_data = False
        # Set up cache
        self.cache_path = tempfile.mkdtemp()
        self.cache = pi.MatrixElementCache(self.cache_path)
        # Setup states
        self.state_one = pi.StateOne("Rb", 61, 2, 1.5, 1.5)
        self.state_two = pi.StateTwo(self.state_one, self.state_one)

    def tearDown(self):
        try:
            shutil.rmtree(self.cache_path)
        except:
            pass

        if self.dump_new_reference_data:
            with open("integration_test_referencedata.pickle", "wb") as f:
                pickle.dump((self.hamiltonian_one, self.basis_one, self.hamiltonian_two, self.basis_two),f)

        self.assertFalse(self.dump_new_reference_data)

    def test_integration(self):
        # Build one-atom system
        system_one = pi.SystemOne(self.state_one.species, self.cache)
        system_one.restrictEnergy(self.state_one.getEnergy() - 40, self.state_one.getEnergy() + 40)
        system_one.restrictN(self.state_one.n - 1, self.state_one.n + 1)
        system_one.restrictL(self.state_one.l - 1, self.state_one.l + 1)
        system_one.setEfield([0, 0, 0.1])
        system_one.setBfield([0, 0, 1])

        # Check for correct dimensions
        self.assertEqual(system_one.getNumVectors(), 64)
        self.assertEqual(system_one.getNumStates(), 64)

        # Compare current results to the reference data (the results have to be
        # compared before diagonalization as the order of the eigenvectors is not
        # fixed)
        hamiltonian_one = system_one.getHamiltonianmatrix()
        basis_one = system_one.getCoefficients()
        # without pruning, max_diff_hamiltonian might be infinity due to
        # division by zero
        hamiltonian_one.data *= (abs(hamiltonian_one).data > 1e-6)
        hamiltonian_one.eliminate_zeros()
        basis_one.data *= (abs(basis_one).data > 1e-6)
        basis_one.eliminate_zeros()

        if self.dump_new_reference_data:
            self.hamiltonian_one = hamiltonian_one.todense()
            self.basis_one = basis_one.todense()
        else:
            with open("integration_test_referencedata.pickle", "rb") as f:
                hamiltonian_one_reference, basis_one_reference, _, _ = pickle.load(f)
                np.testing.assert_allclose(np.asarray(hamiltonian_one.todense()),
                                           np.asarray(hamiltonian_one_reference),
                                           rtol=1e-6)
                np.testing.assert_allclose(np.asarray(basis_one.todense()),
                                           np.asarray(basis_one_reference),
                                           rtol=1e-9)

        # Diagonalize one-atom system
        system_one.diagonalize()

        # Build one-atom system (for this test, system_one has to be diagonal by
        # itself because diagonalization can lead to different order of
        # eigenvectors)
        system_one = pi.SystemOne(self.state_one.species, self.cache)
        system_one.restrictEnergy(self.state_one.getEnergy() - 40, self.state_one.getEnergy() + 40)
        system_one.restrictN(self.state_one.n - 1, self.state_one.n + 1)
        system_one.restrictL(self.state_one.l - 1, self.state_one.l + 1)

        # Build two-atom system
        system_two = pi.SystemTwo(system_one, system_one, self.cache)
        system_two.restrictEnergy(self.state_two.getEnergy() - 2, self.state_two.getEnergy() + 2)
        system_two.setConservedParityUnderPermutation(pi.ODD)
        system_two.setDistance(6)
        system_two.setAngle(0.9)

        # Check for correct dimensions
        self.assertEqual(system_two.getNumVectors(), 239)
        self.assertEqual(system_two.getNumStates(), 468)

        # Compare current results to the reference data (the results have to be
        # compared before diagonalization as the order of the eigenvectors is not
        # fixed)
        hamiltonian_two = system_two.getHamiltonianmatrix()
        basis_two = system_two.getCoefficients()
        # without pruning, max_diff_hamiltonian might be infinity due to
        # division by zero
        hamiltonian_two.data *= (abs(hamiltonian_two).data > 1e-6)
        hamiltonian_two.eliminate_zeros()
        basis_two.data *= (abs(basis_two).data > 1e-6)
        basis_two.eliminate_zeros()

        if self.dump_new_reference_data:
            self.hamiltonian_two = hamiltonian_two.todense()
            self.basis_two = basis_two.todense()
        else:
            with open("integration_test_referencedata.pickle", "rb") as f:
                _, _, hamiltonian_two_reference, basis_two_reference = pickle.load(f)
                np.testing.assert_allclose(np.asarray(hamiltonian_two.todense()),
                                           np.asarray(hamiltonian_two_reference),
                                           rtol=1e-6)
                np.testing.assert_allclose(np.asarray(basis_two.todense()),
                                           np.asarray(basis_two_reference),
                                           rtol=1e-9)

        # Diagonalize two-atom system
        system_two.diagonalize()

if __name__ == '__main__':
    unittest.main()
