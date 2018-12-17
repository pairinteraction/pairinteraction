import numpy as np
import pickle
import shutil
import tempfile
import unittest

from @LIBNAME@ import picomplex as pi


class RotationTest(unittest.TestCase):

    def setUp(self):
        # Set up rotation angles
        self.alpha = 2.67
        self.beta = 1.32
        self.gamma = 0.83

        # Get geometric interpretation of the rotation angles

        # The coordinate system is rotated by R_z(alpha)*R_y(beta)*R_z(gamma) where e.g. R_z is the elementary rotation
        # R_z(alpha) = [[cos(alpha), -sin(alpha), 0], [sin(alpha), cos(alpha), 0], [0,0,1]]. In the following, we calculate
        # the coordinate representations of the rotated z-axis and y-axis in the unrotated system.

        s, c = np.sin(self.gamma), np.cos(self.gamma)
        self.rotator = np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])
        s, c = np.sin(self.beta), np.cos(self.beta)
        self.rotator = np.dot(np.array([[c, 0, s], [0, 1, 0], [-s, 0, c]]), self.rotator)
        s, c = np.sin(self.alpha), np.cos(self.alpha)
        self.rotator = np.dot(np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]]), self.rotator)
        self.to_z_axis = np.dot(self.rotator, [0, 0, 1])
        self.to_y_axis = np.dot(self.rotator, [0, 1, 0])

        # Set up cache
        self.cache = pi.MatrixElementCache()

        # Setup states
        self.state_one = pi.StateOne("Rb", 61, 2, 1.5, 1.5)
        self.state_two = pi.StateTwo(self.state_one, self.state_one)

        # Build one-atom system
        self.system_one = pi.SystemOne(self.state_one.getSpecies(), self.cache)
        self.system_one.restrictEnergy(self.state_one.getEnergy() - 30, self.state_one.getEnergy() + 30)
        self.system_one.restrictN(self.state_one.getN() - 1, self.state_one.getN() + 1)
        self.system_one.restrictL(self.state_one.getL() - 1, self.state_one.getL() + 1)

    def test_rotation_WignerD(self):
        wignerd = pi.WignerD()
        self.assertAlmostEqual(wignerd(0.5, 0.5, 0.5, 0.21, 3.23, 1.23), -0.0332218 + 0.0291377j, places=6)
        self.assertAlmostEqual(wignerd(0.5, 0.5, -0.5, 0.21, 3.23, 1.23), -0.871892 - 0.4877j, places=6)
        self.assertAlmostEqual(wignerd(1.5, 1.5, 1.5, 0.21, 3.23, 1.23), 0.0000479502 + 0.0000717385j, places=6)
        self.assertAlmostEqual(wignerd(1.5, 1.5, 0.5, 0.21, 3.23, 1.23), -0.00202 + 0.00270856j, places=6)
        self.assertAlmostEqual(wignerd(0.5, 0.5, 0.5, np.pi / 2), np.cos(np.pi / 2 / 2))
        self.assertAlmostEqual(wignerd(0.5, 0.5, -0.5, np.pi / 2), -np.sin(np.pi / 2 / 2))
        self.assertAlmostEqual(wignerd(1.5, 1.5, 1.5, np.pi / 2), (1 + np.cos(np.pi / 2)) / 2 * np.cos(np.pi / 2 / 2))
        self.assertAlmostEqual(wignerd(1.5, 1.5, 0.5, np.pi / 2), -np.sqrt(3) *
                               (1 + np.cos(np.pi / 2)) / 2 * np.sin(np.pi / 2 / 2))

    def test_rotation_derotate(self):
        # Add interaction to the Hamiltonian and diagonalize it
        system_one_interacting = pi.SystemOne(self.system_one)
        system_one_interacting.setEfield([0, 0, 0.1])
        system_one_interacting.setBfield([0, 1, 0])
        system_one_interacting.diagonalize()

        # Original eigen pairs
        eigenvectors1 = system_one_interacting.getBasisvectors()
        eigenvalues1 = system_one_interacting.getHamiltonian().diagonal()

        # Eigen pairs after rotating and derotating
        system_one_tmp = pi.SystemOne(system_one_interacting)
        system_one_tmp.rotate(self.alpha, self.beta, self.gamma)
        system_one_tmp.rotate(-self.gamma, -self.beta, -self.alpha)
        eigenvectors2 = system_one_tmp.getBasisvectors()
        eigenvalues2 = system_one_tmp.getHamiltonian().diagonal()

        # Check that everything is the same
        np.testing.assert_allclose(eigenvalues1, eigenvalues2, rtol=1e-6)
        np.testing.assert_allclose(eigenvectors1.A, eigenvectors2.A, rtol=1e-6)

    def test_rotation_overlap(self):
        states_to_calculate_overlap_with = [
            pi.StateOne("Rb", 61, 2, pi.ARB, 1.5),
            pi.StateOne("Rb", 60, 2, 1.5, 1.5)
        ]

        # Add interaction to the Hamiltonian and diagonalize it
        system_one_interacting = pi.SystemOne(self.system_one)
        system_one_interacting.setEfield([0, 0, 0.1])
        system_one_interacting.setBfield([0, 1, 0])
        system_one_interacting.diagonalize()

        # Overlap with rotated system
        system_one_tmp = pi.SystemOne(system_one_interacting)
        system_one_tmp.rotate(self.to_z_axis, self.to_y_axis)
        overlap_rotated_system1 = system_one_tmp.getOverlap(states_to_calculate_overlap_with)

        system_one_tmp = pi.SystemOne(system_one_interacting)
        system_one_tmp.rotate(self.alpha, self.beta, self.gamma)
        overlap_rotated_system2 = system_one_tmp.getOverlap(states_to_calculate_overlap_with)

        system_one_tmp = pi.SystemOne(system_one_interacting)
        # gamma can be arbitrary, it just leads to a phase
        system_one_tmp.rotate(self.alpha, self.beta, self.gamma + 0.8342)
        overlap_rotated_system3 = system_one_tmp.getOverlap(states_to_calculate_overlap_with)

        # Overlap with rotated states
        system_one_tmp = pi.SystemOne(system_one_interacting)
        overlap_rotated_states1 = system_one_tmp.getOverlap(
            states_to_calculate_overlap_with, -self.gamma, -self.beta, -self.alpha)

        system_one_tmp = pi.SystemOne(system_one_interacting)
        overlap_rotated_states2 = system_one_tmp.getOverlap(
            states_to_calculate_overlap_with, -self.gamma + 0.8342, -self.beta, -self.alpha)

        # Check that everything is the same
        np.testing.assert_allclose(overlap_rotated_system1, overlap_rotated_system2, rtol=1e-6)
        np.testing.assert_allclose(overlap_rotated_system1, overlap_rotated_system3, rtol=1e-6)
        np.testing.assert_allclose(overlap_rotated_system1, overlap_rotated_states1, rtol=1e-6)
        np.testing.assert_allclose(overlap_rotated_system1, overlap_rotated_states2, rtol=1e-6)

    def test_rotation_phase(self):
        angle = 0.84

        # Phase due to rotation
        phases1 = self.system_one.getBasisvectors().todense()

        system_one_tmp = pi.SystemOne(self.system_one)
        system_one_tmp.rotate(angle, 0, 0)
        phases2 = system_one_tmp.getBasisvectors().todense()

        # Check diagonality
        np.testing.assert_allclose(phases1, np.diag(np.diag(phases1)), rtol=1e-4, atol=1e-6)
        np.testing.assert_allclose(phases2, np.diag(np.diag(phases2)), rtol=1e-4, atol=1e-6)

        # Check absolute value
        np.testing.assert_allclose(np.abs(phases1), np.abs(phases2), rtol=1e-4, atol=1e-6)

        # Check phases
        array_M = (np.angle(np.diag(phases2)) - np.angle(np.diag(phases1))) / angle
        for state, M in zip(system_one_tmp.getStates(), array_M):
            self.assertAlmostEqual(state.getM(), M)

    def test_rotation_hamiltonian(self):

        # Hamiltonian (in canonical basis) and overlaps after rotating the system with/without atom-field interactions
        system_one_tmp = pi.SystemOne(self.system_one)
        system_one_tmp.setEfield([0.2, 0, 0.3])
        system_one_tmp.setBfield([0, 1, 0])
        system_one_tmp.rotate(self.alpha, self.beta, self.gamma)
        hamiltonian1 = system_one_tmp.getHamiltonian().todense()
        system_one_tmp.diagonalize()
        overlaps1 = system_one_tmp.getOverlap(pi.StateOne("Rb", 61, 2, 1.5, 1.5))

        system_one_tmp = pi.SystemOne(self.system_one)
        # is of no importance for the overlaps, but for the hamiltonians
        system_one_tmp.rotate(self.alpha, self.beta, self.gamma)
        system_one_tmp.setEfield(np.dot(self.rotator.T, [0.2, 0, 0.3]))
        system_one_tmp.setBfield(np.dot(self.rotator.T, [0, 1, 0]))
        hamiltonian2 = system_one_tmp.getHamiltonian().todense()
        system_one_tmp.diagonalize()
        overlaps2 = system_one_tmp.getOverlap(pi.StateOne("Rb", 61, 2, 1.5, 1.5))

        system_one_tmp = pi.SystemOne(self.system_one)
        # is of no importance for the overlaps, but for the hamiltonians
        system_one_tmp.rotate(self.alpha, self.beta, self.gamma)
        system_one_tmp.setEfield([0.2, 0, 0.3], self.alpha, self.beta, self.gamma)
        system_one_tmp.setBfield([0, 1, 0], self.alpha, self.beta, self.gamma)
        hamiltonian3 = system_one_tmp.getHamiltonian().todense()
        system_one_tmp.diagonalize()
        overlaps3 = system_one_tmp.getOverlap(pi.StateOne("Rb", 61, 2, 1.5, 1.5))

        system_one_tmp = pi.SystemOne(self.system_one)
        system_one_tmp.setEfield([0.2, 0, 0.3])
        system_one_tmp.setBfield([0, 1, 0])
        system_one_tmp.diagonalize()
        overlaps4 = system_one_tmp.getOverlap(pi.StateOne("Rb", 61, 2, 1.5, 1.5), -self.gamma, -self.beta, -self.alpha)

        # Check that everything is the same
        np.testing.assert_allclose(overlaps1, overlaps2, rtol=1e-4, atol=1e-6)
        np.testing.assert_allclose(overlaps1, overlaps3, rtol=1e-4, atol=1e-6)
        np.testing.assert_allclose(overlaps1, overlaps4, rtol=1e-4, atol=1e-6)
        np.testing.assert_allclose(hamiltonian1, hamiltonian2, rtol=1e-4, atol=1e-6)
        np.testing.assert_allclose(hamiltonian1, hamiltonian3, rtol=1e-4, atol=1e-6)

    def test_rotation_dipolar(self):
        theta = np.pi / 3

        # Build two-atom system
        system_one_tmp = pi.SystemOne(self.system_one)
        system_one_tmp.setBfield([0, 0, 1])
        system_one_tmp.diagonalize()
        system_two = pi.SystemTwo(system_one_tmp, system_one_tmp, self.cache)
        system_two.restrictEnergy(self.state_two.getEnergy() - 2, self.state_two.getEnergy() + 2)
        system_two.setConservedParityUnderPermutation(pi.ODD)
        system_two.setDistance(6)
        system_two.setAngle(theta)
        system_two.diagonalize()
        eigenvalues1 = system_two.getHamiltonian().diagonal()
        overlaps1 = system_two.getOverlap(self.state_two)

        system_one_tmp = pi.SystemOne(self.system_one)
        system_one_tmp.setBfield([0, 0, 1], 0, theta, 0)
        system_one_tmp.diagonalize()
        system_two = pi.SystemTwo(system_one_tmp, system_one_tmp, self.cache)
        system_two.restrictEnergy(self.state_two.getEnergy() - 2, self.state_two.getEnergy() + 2)
        system_two.setConservedParityUnderPermutation(pi.ODD)
        system_two.setDistance(6)
        system_two.diagonalize()
        eigenvalues2 = system_two.getHamiltonian().diagonal()
        overlaps2 = system_two.getOverlap(self.state_two, 0, theta, 0)

        system_two.rotate(0, -theta, 0)
        overlaps3 = system_two.getOverlap(self.state_two)

        # Check that everything is the same
        np.testing.assert_allclose(eigenvalues1, eigenvalues2, rtol=1e-6)
        np.testing.assert_allclose(overlaps1, overlaps2, rtol=1e-4, atol=1e-6)
        np.testing.assert_allclose(overlaps1, overlaps3, rtol=1e-4, atol=1e-6)


if __name__ == '__main__':
    unittest.main()
