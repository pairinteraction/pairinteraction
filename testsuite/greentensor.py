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
        np.testing.assert_allclose(hamiltonian_standard.A, hamiltonian_greentensor.A, rtol=1e-6)

    def test_greentensor_surface(self):
        theta = np.pi/2
        interatomic_distance = 10
        distance_to_surface = np.array([2.65/6, 5.29/6, 7.9/6])*interatomic_distance # center of mass distance
        state_one1 = pi.StateOne("Rb", 69, 0, 0.5, 0.5)
        state_one2 = pi.StateOne("Rb", 72, 0, 0.5, 0.5)
    
        # Set up pair state
        state_two = pi.StateTwo(state_one1, state_one2)

        # Set up one-atom system
        system_one = pi.SystemOne(state_one1.getSpecies(), self.cache)
        system_one.restrictEnergy(min(state_one1.getEnergy(),state_one2.getEnergy()) - 30, max(state_one1.getEnergy(),state_one2.getEnergy()) + 30)
        system_one.restrictN(min(state_one1.getN(),state_one2.getN()) - 2, max(state_one1.getN(),state_one2.getN()) + 2)
        system_one.restrictL(min(state_one1.getL(),state_one2.getL()) - 1, max(state_one1.getL(),state_one2.getL()) + 1)
        
        # Set up two-atom system
        system_two = pi.SystemTwo(system_one, system_one, self.cache)
        system_two.restrictEnergy(state_two.getEnergy() - 3, state_two.getEnergy() + 3)

        system_two.setAngle(theta)
        system_two.setDistance(interatomic_distance)
        system_two.enableGreenTensor(True)
        
        # Calculate dispersion coefficients
        system_two.diagonalize()
        idx = np.argmax(system_two.getOverlap(state_two, 0, -theta, 0))
        C6_freespace = (system_two.getHamiltonian().diagonal()[idx]-state_two.getEnergy())*interatomic_distance**6
        
        system_two.setSurfaceDistance(distance_to_surface[0])
        system_two.diagonalize()
        idx = np.argmax(system_two.getOverlap(state_two, 0, -theta, 0))
        C6_1 = (system_two.getHamiltonian().diagonal()[idx]-state_two.getEnergy())*interatomic_distance**6
        
        system_two.setSurfaceDistance(distance_to_surface[1])
        system_two.diagonalize()
        idx = np.argmax(system_two.getOverlap(state_two, 0, -theta, 0))
        C6_2 = (system_two.getHamiltonian().diagonal()[idx]-state_two.getEnergy())*interatomic_distance**6
        
        system_two.setSurfaceDistance(distance_to_surface[2])
        system_two.diagonalize()
        idx = np.argmax(system_two.getOverlap(state_two, 0, -theta, 0))
        C6_3 = (system_two.getHamiltonian().diagonal()[idx]-state_two.getEnergy())*interatomic_distance**6

        # Compare the results against literature
        np.testing.assert_allclose(C6_freespace, -670, atol=20)
        np.testing.assert_allclose(C6_1, -544, atol=20)
        np.testing.assert_allclose(C6_2, -628, atol=20)
        np.testing.assert_allclose(C6_3, -649, atol=20)

if __name__ == '__main__':
    unittest.main()
