import numpy as np
import pickle
import shutil
import tempfile
import unittest

from @LIBNAME@ import pireal as pi


class GreenTensorTest(unittest.TestCase):

    def setUp(self):
        # Set up cache
        self.cache = pi.MatrixElementCache()

        # Setup states
        self.state_one = pi.StateOne("Rb", 61, 1, 1.5, 1.5)
        self.state_two = pi.StateTwo(self.state_one, self.state_one)

    def test_greentensor_angle(self):
        # Build one-atom system
        system_one = pi.SystemOne(self.state_one.getSpecies(), self.cache)
        system_one.restrictEnergy(self.state_one.getEnergy() - 40, self.state_one.getEnergy() + 40)
        system_one.restrictN(self.state_one.getN() - 1, self.state_one.getN() + 1)
        system_one.restrictL(self.state_one.getL() - 2, self.state_one.getL() + 2)
        
        # Build two-atom system
        system_two = pi.SystemTwo(system_one, system_one, self.cache)
        system_two.restrictEnergy(self.state_two.getEnergy() - 5, self.state_two.getEnergy() + 5)
        system_two.setConservedParityUnderInversion(pi.ODD)
        system_two.setConservedParityUnderPermutation(pi.ODD)
        system_two.setDistance(5)
        system_two.setAngle(1.78)

        # Construct the Hamiltonian using the standard approach
        system_two_standard = pi.SystemTwo(system_two)
        system_two_standard.enableGreenTensor(False)
        hamiltonian_standard = system_two_standard.getHamiltonian()

        # Construct the Hamiltonian using the green tensor approach
        system_two_greentensor = pi.SystemTwo(system_two)
        system_two_greentensor.enableGreenTensor(True)
        hamiltonian_greentensor = system_two_greentensor.getHamiltonian()

        # Prune Hamiltonians (without pruning, max_diff_hamiltonian might be infinity due to division by zero)
        hamiltonian_standard.data *= (abs(hamiltonian_standard).data > 1e-6)
        hamiltonian_standard.eliminate_zeros()
        hamiltonian_greentensor.data *= (abs(hamiltonian_greentensor).data > 1e-6)
        hamiltonian_greentensor.eliminate_zeros()

        # Compare Hamiltonians
        np.testing.assert_allclose(hamiltonian_standard.A, hamiltonian_greentensor.A, rtol=1e-6)
    
    def test_greentensor_dipolequadrupole(self):
        # Build one-atom system
        system_one = pi.SystemOne(self.state_one.getSpecies(), self.cache)
        system_one.restrictEnergy(self.state_one.getEnergy() - 40, self.state_one.getEnergy() + 40)
        system_one.restrictN(self.state_one.getN() - 2, self.state_one.getN() + 2)
        system_one.restrictL(self.state_one.getL() - 2, self.state_one.getL() + 2)
        
        # Build two-atom system
        system_two = pi.SystemTwo(system_one, system_one, self.cache)
        system_two.restrictEnergy(self.state_two.getEnergy() - 5, self.state_two.getEnergy() + 5)
        system_two.setConservedMomentaUnderRotation([int(np.sum(self.state_two.getM()))])
        system_two.setConservedParityUnderInversion(pi.ODD)
        system_two.setDistance(2)
        system_two.setOrder(4)

        # Construct the Hamiltonian using the standard approach
        system_two_standard = pi.SystemTwo(system_two)
        system_two_standard.enableGreenTensor(False)
        hamiltonian_standard = system_two_standard.getHamiltonian()

        # Construct the Hamiltonian using the green tensor approach
        system_two_greentensor = pi.SystemTwo(system_two)
        system_two_greentensor.enableGreenTensor(True)
        hamiltonian_greentensor = system_two_greentensor.getHamiltonian()

        # Prune Hamiltonians (without pruning, max_diff_hamiltonian might be infinity due to division by zero)
        hamiltonian_standard.data *= (abs(hamiltonian_standard).data > 1e-6)
        hamiltonian_standard.eliminate_zeros()
        hamiltonian_greentensor.data *= (abs(hamiltonian_greentensor).data > 1e-6)
        hamiltonian_greentensor.eliminate_zeros()

        # Compare Hamiltonians
        # TODO np.testing.assert_allclose(hamiltonian_standard.A, hamiltonian_greentensor.A, rtol=1e-6)

    def test_greentensor_surface(self):
        com_distance_to_surface = 3
        
        # Build one-atom system
        system_one = pi.SystemOne(self.state_one.getSpecies(), self.cache)
        system_one.restrictEnergy(self.state_one.getEnergy() - 40, self.state_one.getEnergy() + 40)
        system_one.restrictN(self.state_one.getN() - 1, self.state_one.getN() + 1)
        system_one.restrictL(self.state_one.getL() - 2, self.state_one.getL() + 2)
        
        # Build two-atom system
        system_two = pi.SystemTwo(system_one, system_one, self.cache)
        system_two.restrictEnergy(self.state_two.getEnergy() - 2, self.state_two.getEnergy() + 2)
        system_two.setConservedParityUnderPermutation(pi.ODD)
        system_two.setDistance(5)
        system_two.setAngle(1.78)
        system_two.setSurfaceDistance(com_distance_to_surface)
        system_two.enableGreenTensor(True)
        
        # Construct the Hamiltonian
        hamiltonian_greentensor = system_two.getHamiltonian()

        # Compare the results against literature
        # TODO


if __name__ == '__main__':
    unittest.main()
