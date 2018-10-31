import numpy as np
import pickle
import shutil
import tempfile
import unittest

from @LIBNAME@ import pireal as pi


class IntegrationTest(unittest.TestCase):

    def setUp(self):
        self.cache_path = tempfile.mkdtemp()
        self.cache = pi.MatrixElementCache(self.cache_path)

    def tearDown(self):
        try:
            shutil.rmtree(self.cache_path)
        except BaseException:
            pass

    def test_perturbation_zero_angle(self):

        distance = 5
        deltaN = 2

        state_two_subspace = [pi.StateTwo(["Rb", "Rb"], [42, 42], [0, 1], [1 / 2, 1 / 2], [1 / 2, -1 / 2]),
                              pi.StateTwo(["Rb", "Rb"], [42, 42], [0, 1], [1 / 2, 1 / 2], [-1 / 2, 1 / 2]),
                              pi.StateTwo(["Rb", "Rb"], [42, 42], [1, 0], [1 / 2, 1 / 2], [1 / 2, -1 / 2]),
                              pi.StateTwo(["Rb", "Rb"], [42, 42], [1, 0], [1 / 2, 1 / 2], [-1 / 2, 1 / 2])]

        ### Schrieffer Wolff transformation ###

        state_one = state_two_subspace[0].getFirstState()
        system_one = pi.SystemOne(state_one.getSpecies(), self.cache)
        system_one.restrictEnergy(state_one.getEnergy() - 200, state_one.getEnergy() + 200)
        system_one.restrictN(state_one.getN() - deltaN, state_one.getN() + deltaN)
        system_one.restrictL(state_one.getL() - 1, state_one.getL() + 1)
        system_two = pi.SystemTwo(system_one, system_one, self.cache)
        system_two.restrictEnergy(state_two_subspace[0].getEnergy() - 20, state_two_subspace[0].getEnergy() + 20)
        system_two.setConservedMomentaUnderRotation([int(np.sum(state_two_subspace[0].getM()))])

        # System containing the unperturbed Hamiltonian
        system_two_unperturbed = pi.SystemTwo(system_two)

        # System containing the perturbed Hamiltonian
        system_two_perturbed = pi.SystemTwo(system_two)
        system_two_perturbed.setDistance(distance)

        # Restrict unperturbed system to the subspace
        system_two_unperturbed.constrainBasisvectors(system_two_unperturbed.getBasisvectorIndex(state_two_subspace))

        # Apply the Schrieffer Wolff transformation on the perturbed Hamiltonian
        system_two_perturbed.applySchriefferWolffTransformation(system_two_unperturbed)

        # Effective Hamiltonian
        hamiltonian_schrieffer_wolff = system_two_perturbed.getHamiltonian()

        ### Perturbation theory up to second order ###

        calculator = pi.PerturbativeInteraction(self.cache)

        E0 = calculator.getEnergy(state_two_subspace)
        C3 = calculator.getC3(state_two_subspace)
        C6 = calculator.getC6(state_two_subspace, deltaN)

        # Effective Hamiltonian
        hamiltonian_second_order = E0 + C3 / distance**3 + C6 / distance**6

        ### Compare results ###

        np.testing.assert_allclose(hamiltonian_schrieffer_wolff.A, hamiltonian_second_order, rtol=1e-2)

    def test_perturbation_nonzero_angle(self):

        theta = np.pi / 3
        distance = 5
        deltaN = 2

        state_two_subspace = [pi.StateTwo(["Rb", "Rb"], [42, 42], [0, 0], [1 / 2, 1 / 2], [-1 / 2, -1 / 2]),
                              pi.StateTwo(["Rb", "Rb"], [42, 42], [0, 0], [1 / 2, 1 / 2], [-1 / 2, 1 / 2]),
                              pi.StateTwo(["Rb", "Rb"], [42, 42], [0, 0], [1 / 2, 1 / 2], [1 / 2, -1 / 2]),
                              pi.StateTwo(["Rb", "Rb"], [42, 42], [0, 0], [1 / 2, 1 / 2], [1 / 2, 1 / 2])]

        ### Schrieffer Wolff transformation ###

        state_one = state_two_subspace[0].getFirstState()
        system_one = pi.SystemOne(state_one.getSpecies(), self.cache)
        system_one.restrictEnergy(state_one.getEnergy() - 200, state_one.getEnergy() + 200)
        system_one.restrictN(state_one.getN() - deltaN, state_one.getN() + deltaN)
        system_one.restrictL(state_one.getL() - 1, state_one.getL() + 1)
        system_two = pi.SystemTwo(system_one, system_one, self.cache)
        system_two.restrictEnergy(state_two_subspace[0].getEnergy() - 20, state_two_subspace[0].getEnergy() + 20)

        # System containing the unperturbed Hamiltonian
        system_two_unperturbed = pi.SystemTwo(system_two)

        # System containing the perturbed Hamiltonian
        system_two_perturbed = pi.SystemTwo(system_two)
        system_two_perturbed.setDistance(distance)
        system_two_perturbed.setAngle(theta)

        # Restrict unperturbed system to the subspace
        system_two_unperturbed.constrainBasisvectors(system_two_unperturbed.getBasisvectorIndex(state_two_subspace))

        # Apply the Schrieffer Wolff transformation on the perturbed Hamiltonian
        system_two_perturbed.applySchriefferWolffTransformation(system_two_unperturbed)

        # Effective Hamiltonian
        hamiltonian_schrieffer_wolff = system_two_perturbed.getHamiltonian()

        ### Perturbation theory up to second order ###

        calculator = pi.PerturbativeInteraction(theta, self.cache)

        E0 = calculator.getEnergy(state_two_subspace)
        C3 = calculator.getC3(state_two_subspace)
        C6 = calculator.getC6(state_two_subspace, deltaN)

        # Effective Hamiltonian
        hamiltonian_second_order = E0 + C3 / distance**3 + C6 / distance**6

        ### Compare results ###

        np.testing.assert_allclose(hamiltonian_schrieffer_wolff.A, hamiltonian_second_order, rtol=1e-2)

    def test_perturbation_C6(self):
        calculator = pi.PerturbativeInteraction(self.cache)
        C6 = calculator.getC6(pi.StateTwo(["Rb", "Rb"], [60, 60], [2, 2], [3 / 2, 3 / 2], [3 / 2, 3 / 2]), 4)

        self.assertAlmostEqual(C6, -886.744, places=3)


if __name__ == '__main__':
    unittest.main()
