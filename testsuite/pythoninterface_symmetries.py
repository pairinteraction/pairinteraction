import unittest
import numpy as np
from @LIBNAME@ import picomplex as pi
import tempfile
import shutil

class TestPythoninterfaceSymmetries(unittest.TestCase):

    def test_comparison(self):

        #######################################################
        ### Preparations ######################################
        #######################################################

        # Create cache directory
        path_cache = tempfile.mkdtemp()

        #######################################################
        ### Check rotation symmetry of one atom systems #######
        #######################################################

        # Define state
        state_one = pi.StateOne("Rb", 61, 1, 0.5, 0.5)

        # Build one atom system
        system_one = pi.SystemOne(state_one.element, path_cache)
        system_one.restrictEnergy(state_one.energy-40, state_one.energy+40)
        system_one.restrictN(state_one.n-1, state_one.n+1)
        system_one.restrictL(state_one.l-1, state_one.l+1)
        system_one.restrictJ(state_one.j-1, state_one.j+1)
        system_one.setBfield([0, 0, 100])
        system_one.setEfield([0, 0, 1])

        momenta_one = np.arange(-(state_one.j+1), (state_one.j+1)+1)

        # Diagonalize blockwise
        system_one_momentum = {}
        for m in momenta_one:
            system_one_momentum[m] = pi.SystemOne(system_one)
            system_one_momentum[m].setConservedMomentaUnderRotation([m])
            system_one_momentum[m].diagonalize()

        system_one_combined = pi.SystemOne(system_one_momentum[momenta_one[0]])
        for m in momenta_one[1:]:
            system_one_combined.add(system_one_momentum[m])

        # Diagonalize altogether
        system_one_alternative = pi.SystemOne(system_one)
        system_one_alternative.setConservedMomentaUnderRotation(momenta_one)
        system_one_alternative.diagonalize()

        system_one.setConservedMomentaUnderRotation([pi.ARB])
        system_one.diagonalize()

        # Compare results
        w1 = np.sort(system_one_combined.diagonal)
        w2 = np.sort(system_one.diagonal)
        w3 = np.sort(system_one_alternative.diagonal)

        maxdiff12 = np.abs((w1-w2)/(np.max([w1,w2],axis=0)))
        maxdiff23 = np.abs((w3-w2)/(np.max([w3,w2],axis=0)))
        maxdiff12[np.abs(w1-w2)<1e-14] = np.abs(w1-w2)[np.abs(w1-w2)<1e-14]
        maxdiff23[np.abs(w3-w2)<1e-14] = np.abs(w3-w2)[np.abs(w3-w2)<1e-14]
        maxdiff12 = np.max(maxdiff12)
        maxdiff23 = np.max(maxdiff23)
        print("One-atom system with rotation symmetry, relative maximum deviation: ", maxdiff12, " (between alternatives: ", maxdiff23,")")
        self.assertAlmostEqual(maxdiff12, 0, places=9)
        self.assertAlmostEqual(maxdiff23, 0, places=9)

        #######################################################
        ### Check rotation symmetry of two atom systems #######
        #######################################################

        # Define state
        state_two = pi.StateTwo(state_one, state_one)

        # Diagonalize blockwise
        system_two_momentum = {}

        for m1 in momenta_one:
            for m2 in momenta_one:
                if m1+m2 in system_two_momentum:
                    tmp = pi.SystemTwo(system_one_momentum[m1], system_one_momentum[m2])
                    tmp.restrictEnergy(state_two.energy-2, state_two.energy+2)
                    system_two_momentum[m1+m2].add(tmp)
                else:
                    system_two_momentum[m1+m2] = pi.SystemTwo(system_one_momentum[m1], system_one_momentum[m2])
                    system_two_momentum[m1+m2].restrictEnergy(state_two.energy-2, state_two.energy+2)

        momenta_two = list(system_two_momentum.keys())

        for m in momenta_two:
            system_two_momentum[m].setDistance(1)
            system_two_momentum[m].setOrder(5)
            system_two_momentum[m].diagonalize()

        system_two_combined = pi.SystemTwo(system_two_momentum[momenta_two[0]])
        for m in momenta_two[1:]:
            system_two_combined.add(system_two_momentum[m])

        # Diagonalize blockwise, alternative
        system_two_momentum_alternative = {}
        for m in momenta_two:
            system_two_momentum_alternative[m] = pi.SystemTwo(system_one, system_one)
            system_two_momentum_alternative[m].restrictEnergy(state_two.energy-2, state_two.energy+2)
            system_two_momentum_alternative[m].setConservedMomentaUnderRotation([int(m)])
            system_two_momentum_alternative[m].setDistance(1)
            system_two_momentum_alternative[m].setOrder(5)
            system_two_momentum_alternative[m].diagonalize()

        system_two_combined_alternative = pi.SystemTwo(system_two_momentum_alternative[momenta_two[0]])
        for m in momenta_two[1:]:
            system_two_combined_alternative.add(system_two_momentum_alternative[m])

        # Diagonalize altogether
        system_two = pi.SystemTwo(system_one, system_one, path_cache)
        system_two.restrictEnergy(state_two.energy-2, state_two.energy+2)
        system_two.setDistance(1)
        system_two.setOrder(5)
        system_two.diagonalize()

        # Compare results
        w1 = np.sort(system_two_combined.diagonal)
        w2 = np.sort(system_two.diagonal)
        w3 = np.sort(system_two_combined_alternative.diagonal)

        maxdiff12 = np.abs((w1-w2)/(np.max([w1,w2],axis=0)))
        maxdiff23 = np.abs((w3-w2)/(np.max([w3,w2],axis=0)))
        maxdiff12[np.abs(w1-w2)<1e-14] = np.abs(w1-w2)[np.abs(w1-w2)<1e-14]
        maxdiff23[np.abs(w3-w2)<1e-14] = np.abs(w3-w2)[np.abs(w3-w2)<1e-14]
        maxdiff12 = np.max(maxdiff12)
        maxdiff23 = np.max(maxdiff23)
        print("Two-atom system with rotation symmetry, relative maximum deviation: ", maxdiff12, " (between alternatives: ", maxdiff23,")")
        self.assertAlmostEqual(maxdiff12, 0, places=9)
        self.assertAlmostEqual(maxdiff23, 0, places=9)

        #######################################################
        ### Check reflection symmetry of one atom systems #####
        #######################################################

        # Define state
        state_one = pi.StateOne("Rb", 61, 1, 0.5, 0.5)

        # Build one atom system
        system_one = pi.SystemOne(state_one.element, path_cache)
        system_one.restrictEnergy(state_one.energy-20, state_one.energy+20)
        system_one.restrictN(state_one.n-1, state_one.n+1)
        system_one.restrictL(state_one.l-1, state_one.l+1)
        system_one.restrictJ(state_one.j-1, state_one.j+1)
        system_one.setBfield([0, 100, 0])
        system_one.setEfield([1, 0, 2])

        # Diagonalize blockwise
        system_one_even = pi.SystemOne(system_one)
        system_one_even.setConservedParityUnderReflection(pi.EVEN)
        system_one_even.diagonalize()

        system_one_odd = pi.SystemOne(system_one)
        system_one_odd.setConservedParityUnderReflection(pi.ODD)
        system_one_odd.diagonalize()

        system_one_combined = pi.SystemOne(system_one_even)
        system_one_combined.add(system_one_odd)

        # Diagonalize altogether
        system_one.setConservedParityUnderReflection(pi.NA)
        system_one.diagonalize()

        # Compare results
        w1 = np.sort(system_one_combined.diagonal)
        w2 = np.sort(system_one.diagonal)

        maxdiff = np.abs((w1-w2)/(np.max([w1,w2],axis=0)))
        maxdiff[np.abs(w1-w2)<1e-14] = np.abs(w1-w2)[np.abs(w1-w2)<1e-14]
        maxdiff = np.max(maxdiff)
        print("One-atom system with reflection symmetry, relative maximum deviation: ", maxdiff)
        self.assertAlmostEqual(maxdiff, 0, places=9)

        #######################################################
        ### Check reflection symmetry of two atom systems #####
        #######################################################

        # Remark: calling restrictEnergy() would cause a small deviation

        # Define state
        state_two = pi.StateTwo(state_one, state_one)

        # Diagonalize blockwise
        system_two_even = pi.SystemTwo(system_one_even, system_one_even, path_cache)
        system_two_even.add(pi.SystemTwo(system_one_odd, system_one_odd, path_cache))
        system_two_even.setDistance(1)
        system_two_even.setOrder(5)
        system_two_even.diagonalize()

        system_two_odd = pi.SystemTwo(system_one_even, system_one_odd, path_cache)
        system_two_odd.add(pi.SystemTwo(system_one_odd, system_one_even, path_cache))
        system_two_odd.setDistance(1)
        system_two_odd.setOrder(5)
        system_two_odd.diagonalize()

        system_two_combined =  pi.SystemTwo(system_two_even)
        system_two_combined.add(system_two_odd)

        # Diagonalize altogether
        system_two = pi.SystemTwo(system_one, system_one, path_cache)
        system_two.setDistance(1)
        system_two.setOrder(5)
        system_two.diagonalize()

        # Compare results
        w1 = np.sort(system_two_combined.diagonal)
        w2 = np.sort(system_two.diagonal)

        maxdiff = np.abs((w1-w2)/(np.max([w1,w2],axis=0)))
        maxdiff[np.abs(w1-w2)<1e-14] = np.abs(w1-w2)[np.abs(w1-w2)<1e-14]
        maxdiff = np.max(maxdiff)
        print("Two-atom system with reflection symmetry, relative maximum deviation: ", maxdiff)
        self.assertAlmostEqual(maxdiff, 0, places=9)

        #######################################################
        ### Check permutation symmetry of two atom systems ####
        #######################################################

        # Remark: calling restrictEnergy() would cause a small deviation # TODO figure out exact reason (in case of symmetrized basis states, the atom-atom interaction would add energy to the diagonal)

        # Define states
        state_one = pi.StateOne("Rb", 61, 1, 0.5, 0.5)
        state_two = pi.StateTwo(state_one, state_one)

        # Build one atom system
        system_one = pi.SystemOne(state_one.element, path_cache)
        system_one.restrictEnergy(state_one.energy-20, state_one.energy+20)
        system_one.restrictN(state_one.n-1, state_one.n+1)
        system_one.restrictL(state_one.l-1, state_one.l+1)
        system_one.restrictJ(state_one.j-1, state_one.j+1)
        system_one.setBfield([100, 200, 300])
        system_one.setEfield([1, 2, 3])

        # Build two atom system
        system_two = pi.SystemTwo(system_one, system_one, path_cache)
        #system_two.restrictEnergy(state_two.energy-2, state_two.energy+2) # TODO
        system_two.setDistance(1)
        system_two.setOrder(3)

        # Diagonalize blockwise
        system_two_even = pi.SystemTwo(system_two)
        system_two_even.setConservedParityUnderPermutation(pi.EVEN)
        system_two_even.diagonalize()

        system_two_odd = pi.SystemTwo(system_two)
        system_two_odd.setConservedParityUnderPermutation(pi.ODD)
        system_two_odd.diagonalize()

        system_two_combined = system_two_even
        system_two_combined.add(system_two_odd)

        # Diagonalize altogether
        system_two.setConservedParityUnderPermutation(pi.NA)
        system_two.diagonalize()

        # Compare results
        w1 = np.sort(system_two_combined.diagonal)
        w2 = np.sort(system_two.diagonal)

        maxdiff = np.abs((w1-w2)/(np.max([w1,w2],axis=0)))
        maxdiff[np.abs(w1-w2)<1e-14] = np.abs(w1-w2)[np.abs(w1-w2)<1e-14]
        maxdiff = np.max(maxdiff)
        print("Two-atom system with permutation symmetry, relative maximum deviation: ", maxdiff)
        self.assertAlmostEqual(maxdiff, 0, places=9)

        #######################################################
        ### Check inversion symmetry of two atom systems ######
        #######################################################

        # Remark: calling restrictEnergy() would cause a small deviation # TODO figure out exact reason (in case of symmetrized basis states, the atom-atom interaction would add energy to the diagonal)

        # Define states
        state_one = pi.StateOne("Rb", 61, 1, 0.5, 0.5)
        state_two = pi.StateTwo(state_one, state_one)

        # Build one atom system
        system_one = pi.SystemOne(state_one.element, path_cache)
        system_one.restrictEnergy(state_one.energy-20, state_one.energy+20)
        system_one.restrictN(state_one.n-1, state_one.n+1)
        system_one.restrictL(state_one.l-1, state_one.l+1)
        system_one.restrictJ(state_one.j-1, state_one.j+1)
        system_one.setBfield([100, 200, 300])

        # Build two atom system
        system_two = pi.SystemTwo(system_one, system_one, path_cache)
        #system_two.restrictEnergy(state_two.energy-2, state_two.energy+2) # TODO
        system_two.setDistance(1)
        system_two.setOrder(5)

        # Diagonalize blockwise
        system_two_even = pi.SystemTwo(system_two)
        system_two_even.setConservedParityUnderInversion(pi.EVEN)
        system_two_even.diagonalize()

        system_two_odd = pi.SystemTwo(system_two)
        system_two_odd.setConservedParityUnderInversion(pi.ODD)
        system_two_odd.diagonalize()

        system_two_combined = system_two_even
        system_two_combined.add(system_two_odd)

        # Diagonalize altogether
        system_two.setConservedParityUnderInversion(pi.NA)
        system_two.diagonalize()

        # Compare results
        w1 = np.sort(system_two_combined.diagonal)
        w2 = np.sort(system_two.diagonal)

        maxdiff = np.abs((w1-w2)/(np.max([w1,w2],axis=0)))
        maxdiff[np.abs(w1-w2)<1e-14] = np.abs(w1-w2)[np.abs(w1-w2)<1e-14]
        maxdiff = np.max(maxdiff)
        print("Two-atom system with inversion symmetry, relative maximum deviation: ", maxdiff)
        self.assertAlmostEqual(maxdiff, 0, places=9)

        """#######################################################
        ### Check combined binary symmetries ##################
        #######################################################

        # Define states
        state_one = pi.StateOne("Rb", 61, 1, 0.5, 0.5)
        state_two = pi.StateTwo(state_one, state_one)

        # Build one atom system
        system_one = pi.SystemOne(state_one.element, path_cache)
        system_one.restrictEnergy(state_one.energy-20, state_one.energy+20)
        system_one.restrictN(state_one.n-1, state_one.n+1)
        system_one.restrictL(state_one.l-1, state_one.l+1)
        system_one.restrictJ(state_one.j-1, state_one.j+1)
        system_one.setBfield([0, 100, 0])

        system_one_even = pi.SystemOne(system_one)
        system_one_even.setConservedParityUnderReflection(pi.EVEN)

        system_one_odd = pi.SystemOne(system_one)
        system_one_odd.setConservedParityUnderReflection(pi.ODD)

        system_one.setConservedParityUnderReflection(pi.NA)

        # Build two atom system
        system_two_even = pi.SystemTwo(system_one_even, system_one_even, path_cache)
        system_two_even.add(pi.SystemTwo(system_one_odd, system_one_odd, path_cache))
        #system_two_even.restrictEnergy(state_two.energy-2, state_two.energy+2) # TODO
        system_two_even.setDistance(1)
        system_two_even.setOrder(5)

        system_two_odd = pi.SystemTwo(system_one_even, system_one_odd, path_cache)
        system_two_odd.add(pi.SystemTwo(system_one_odd, system_one_even, path_cache))
        #system_two_odd.restrictEnergy(state_two.energy-2, state_two.energy+2) # TODO
        system_two_odd.setDistance(1)
        system_two_odd.setOrder(5)

        system_two = pi.SystemTwo(system_one, system_one, path_cache)
        #system_two.restrictEnergy(state_two.energy-2, state_two.energy+2) # TODO
        system_two.setDistance(1)
        system_two.setOrder(5)

        # Diagonalize blockwise
        system_two_even_even_even = pi.SystemTwo(system_two_even)
        system_two_even_even_even.setConservedParityUnderPermutation(pi.EVEN)
        system_two_even_even_even.setConservedParityUnderInversion(pi.EVEN)
        system_two_even_even_even.diagonalize()

        system_two_even_even_odd = pi.SystemTwo(system_two_even)
        system_two_even_even_odd.setConservedParityUnderPermutation(pi.EVEN)
        system_two_even_even_odd.setConservedParityUnderInversion(pi.ODD)
        system_two_even_even_odd.diagonalize()

        system_two_even_odd_even = pi.SystemTwo(system_two_even)
        system_two_even_odd_even.setConservedParityUnderPermutation(pi.ODD)
        system_two_even_odd_even.setConservedParityUnderInversion(pi.EVEN)
        system_two_even_odd_even.diagonalize()

        system_two_even_odd_odd = pi.SystemTwo(system_two_even)
        system_two_even_odd_odd.setConservedParityUnderPermutation(pi.ODD)
        system_two_even_odd_odd.setConservedParityUnderInversion(pi.ODD)
        system_two_even_odd_odd.diagonalize()

        system_two_odd_even_even = pi.SystemTwo(system_two_odd)
        system_two_odd_even_even.setConservedParityUnderPermutation(pi.EVEN)
        system_two_odd_even_even.setConservedParityUnderInversion(pi.EVEN)
        system_two_odd_even_even.diagonalize()

        system_two_odd_even_odd = pi.SystemTwo(system_two_odd)
        system_two_odd_even_odd.setConservedParityUnderPermutation(pi.EVEN)
        system_two_odd_even_odd.setConservedParityUnderInversion(pi.ODD)
        system_two_odd_even_odd.diagonalize()

        system_two_odd_odd_even = pi.SystemTwo(system_two_odd)
        system_two_odd_odd_even.setConservedParityUnderPermutation(pi.ODD)
        system_two_odd_odd_even.setConservedParityUnderInversion(pi.EVEN)
        system_two_odd_odd_even.diagonalize()

        system_two_odd_odd_odd = pi.SystemTwo(system_two_odd)
        system_two_odd_odd_odd.setConservedParityUnderPermutation(pi.ODD)
        system_two_odd_odd_odd.setConservedParityUnderInversion(pi.ODD)
        system_two_odd_odd_odd.diagonalize()

        system_two_combined = system_two_even_even_even
        system_two_combined.add(system_two_even_even_odd)
        system_two_combined.add(system_two_even_odd_even)
        system_two_combined.add(system_two_even_odd_odd)
        system_two_combined.add(system_two_odd_even_even)
        system_two_combined.add(system_two_odd_even_odd)
        system_two_combined.add(system_two_odd_odd_even)
        system_two_combined.add(system_two_odd_odd_odd)

        # Diagonalize altogether
        system_two.setConservedParityUnderPermutation(pi.NA)
        system_two.setConservedParityUnderInversion(pi.NA)
        system_two.diagonalize()

        # Compare results
        w1 = np.sort(system_two_combined.diagonal)
        w2 = np.sort(system_two.diagonal)

        maxdiff = np.abs((w1-w2)/(np.max([w1,w2],axis=0)))
        maxdiff[np.abs(w1-w2)<1e-14] = np.abs(w1-w2)[np.abs(w1-w2)<1e-14]
        maxdiff = np.max(maxdiff)
        print("Two-atom system with combined binary symmetries, relative maximum deviation: ", maxdiff)
        self.assertAlmostEqual(maxdiff, 0, places=9)""" # TODO make this work, circumventing error "One cannot change symmetries after the basis was built." (system_two_odd.add builds the basis)

        #######################################################
        ### Clean up ##########################################
        #######################################################

        # Delete cache directory
        shutil.rmtree(path_cache)

if __name__ == '__main__':
    unittest.main()
