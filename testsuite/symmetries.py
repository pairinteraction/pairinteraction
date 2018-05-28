import unittest
import numpy as np
from @LIBNAME@ import picomplex as pi
import tempfile
import shutil
from itertools import product

class TestPythoninterfaceSymmetries(unittest.TestCase):

    #######################################################
    ### Preparations ######################################
    #######################################################
    def setUp(self):
        # Create cache directory
        self.path_cache = tempfile.mkdtemp()

        # Set up cache
        self.cache = pi.MatrixElementCache(self.path_cache)

    def test_rotation(self):
        #######################################################
        ### Check rotation symmetry of one atom systems #######
        #######################################################

        # Define state
        state_one = pi.StateOne("Rb", 61, 1, 0.5, 0.5)

        # Build one atom system
        system_one = pi.SystemOne(state_one.species, self.cache)
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
        print("One-atom system with rotation symmetry, relative maximum deviation: ",
              maxdiff12, " (between alternatives: ", maxdiff23,")")
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
                    tmp = pi.SystemTwo(system_one_momentum[m1], system_one_momentum[m2], self.cache)
                    tmp.restrictEnergy(state_two.energy-2, state_two.energy+2)
                    system_two_momentum[m1+m2].add(tmp)
                else:
                    system_two_momentum[m1+m2] = pi.SystemTwo(system_one_momentum[m1], system_one_momentum[m2], self.cache)
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
            system_two_momentum_alternative[m] = pi.SystemTwo(system_one, system_one, self.cache)
            system_two_momentum_alternative[m].restrictEnergy(state_two.energy-2, state_two.energy+2)
            system_two_momentum_alternative[m].setConservedMomentaUnderRotation([int(m)])
            system_two_momentum_alternative[m].setDistance(1)
            system_two_momentum_alternative[m].setOrder(5)
            system_two_momentum_alternative[m].diagonalize()

        system_two_combined_alternative = pi.SystemTwo(system_two_momentum_alternative[momenta_two[0]])
        for m in momenta_two[1:]:
            system_two_combined_alternative.add(system_two_momentum_alternative[m])

        # Diagonalize altogether
        system_two = pi.SystemTwo(system_one, system_one, self.cache)
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
        print("Two-atom system with rotation symmetry, relative maximum deviation: ",
              maxdiff12, " (between alternatives: ", maxdiff23,")")
        self.assertAlmostEqual(maxdiff12, 0, places=9)
        self.assertAlmostEqual(maxdiff23, 0, places=9)

    def test_reflection(self):
        #######################################################
        ### Check reflection symmetry of one atom systems #####
        #######################################################

        # Define state
        state_one = pi.StateOne("Rb", 61, 1, 0.5, 0.5)

        # Build one atom system
        system_one = pi.SystemOne(state_one.species, self.cache)
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
        # Note: it is called odd in order to fit to the notion of the paper
        system_two_odd = pi.SystemTwo(system_one_even, system_one_even, self.cache)
        system_two_odd.add(pi.SystemTwo(system_one_odd, system_one_odd, self.cache))
        system_two_odd.setDistance(1)
        system_two_odd.setOrder(5)
        system_two_odd.diagonalize()

        system_two_even = pi.SystemTwo(system_one_even, system_one_odd, self.cache)
        system_two_even.add(pi.SystemTwo(system_one_odd, system_one_even, self.cache))
        system_two_even.setDistance(1)
        system_two_even.setOrder(5)
        system_two_even.diagonalize()

        system_two_combined =  pi.SystemTwo(system_two_even)
        system_two_combined.add(system_two_odd)

        # Diagonalize blockwise alternative
        # Note: it is important to use system_one_combined
        system_two_alternative= pi.SystemTwo(system_one_combined, system_one_combined, self.cache)
        system_two_alternative.setDistance(1)
        system_two_alternative.setOrder(5)

        system_two_even_alternative = pi.SystemTwo(system_two_alternative)
        system_two_even_alternative.setConservedParityUnderReflection(pi.EVEN)
        system_two_even_alternative.diagonalize()

        system_two_odd_alternative = pi.SystemTwo(system_two_alternative)
        system_two_odd_alternative.setConservedParityUnderReflection(pi.ODD)
        system_two_odd_alternative.diagonalize()

        # Diagonalize altogether
        system_two = pi.SystemTwo(system_one, system_one, self.cache)
        system_two.setDistance(1)
        system_two.setOrder(5)
        system_two.diagonalize()

        # Compare results
        w1 = np.sort(system_two_combined.diagonal)
        w2 = np.sort(system_two.diagonal)
        w3 = np.sort(system_two_even.diagonal)
        w4 = np.sort(system_two_odd.diagonal)
        w5 = np.sort(system_two_even_alternative.diagonal)
        w6 = np.sort(system_two_odd_alternative.diagonal)

        maxdiff12 = np.abs((w1-w2)/(np.max([w1,w2],axis=0)))
        maxdiff35 = np.abs((w3-w5)/(np.max([w3,w5],axis=0)))
        maxdiff46 = np.abs((w4-w6)/(np.max([w4,w6],axis=0)))
        maxdiff12[np.abs(w1-w2)<1e-14] = np.abs(w1-w2)[np.abs(w1-w2)<1e-14]
        maxdiff35[np.abs(w3-w5)<1e-14] = np.abs(w3-w5)[np.abs(w3-w5)<1e-14]
        maxdiff46[np.abs(w4-w6)<1e-14] = np.abs(w4-w6)[np.abs(w4-w6)<1e-14]
        maxdiff12 = np.max(maxdiff12)
        maxdiff35 = np.max(maxdiff35)
        maxdiff46 = np.max(maxdiff46)
        print("Two-atom system with reflection symmetry, relative maximum deviation: ",
              maxdiff12, " (between alternatives: ", maxdiff35,", ", maxdiff46, ")")
        self.assertAlmostEqual(maxdiff12, 0, places=9)
        self.assertAlmostEqual(maxdiff35, 0, places=9)
        self.assertAlmostEqual(maxdiff46, 0, places=9)

    def test_permutation(self):
        #######################################################
        ### Check permutation symmetry of two atom systems ####
        #######################################################

        # Remark: calling restrictEnergy() would cause a small
        # deviation # TODO figure out exact reason (in case of
        # symmetrized basis states, the atom-atom interaction would
        # add energy to the diagonal)

        # Define states
        state_one = pi.StateOne("Rb", 61, 1, 0.5, 0.5)
        state_two = pi.StateTwo(state_one, state_one)

        # Build one atom system
        system_one = pi.SystemOne(state_one.species, self.cache)
        system_one.restrictEnergy(state_one.energy-20, state_one.energy+20)
        system_one.restrictN(state_one.n-1, state_one.n+1)
        system_one.restrictL(state_one.l-1, state_one.l+1)
        system_one.restrictJ(state_one.j-1, state_one.j+1)
        system_one.setBfield([100, 200, 300])
        system_one.setEfield([1, 2, 3])

        # Build two atom system
        system_two = pi.SystemTwo(system_one, system_one, self.cache)
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

    def test_inversion(self):
        #######################################################
        ### Check inversion symmetry of two atom systems ######
        #######################################################

        # Remark: calling restrictEnergy() would cause a small
        # deviation # TODO figure out exact reason (in case of
        # symmetrized basis states, the atom-atom interaction would
        # add energy to the diagonal)

        # Define states
        state_one = pi.StateOne("Rb", 61, 1, 0.5, 0.5)
        state_two = pi.StateTwo(state_one, state_one)

        # Build one atom system
        system_one = pi.SystemOne(state_one.species, self.cache)
        system_one.restrictEnergy(state_one.energy-20, state_one.energy+20)
        system_one.restrictN(state_one.n-1, state_one.n+1)
        system_one.restrictL(state_one.l-1, state_one.l+1)
        system_one.restrictJ(state_one.j-1, state_one.j+1)
        system_one.setBfield([100, 200, 300])

        # Build two atom system
        system_two = pi.SystemTwo(system_one, system_one, self.cache)
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

    def test_combined(self):
        #######################################################
        ### Check combined binary symmetries ##################
        #######################################################

        # Remark: calling restrictEnergy() would cause a small
        # deviation # TODO figure out exact reason (in case of
        # symmetrized basis states, the atom-atom interaction would
        # add energy to the diagonal)

        # Define states
        state_one = pi.StateOne("Rb", 61, 1, 0.5, 0.5)
        state_two = pi.StateTwo(state_one, state_one)

        # Build one atom system
        system_one = pi.SystemOne(state_one.species, self.cache)
        system_one.restrictEnergy(state_one.energy-20, state_one.energy+20)
        system_one.restrictN(state_one.n-1, state_one.n+1)
        system_one.restrictL(state_one.l-1, state_one.l+1)
        system_one.restrictJ(state_one.j-1, state_one.j+1)
        system_one.setBfield([0, 100, 0])

        system_one_combined = pi.SystemOne(system_one)
        system_one_combined.setConservedParityUnderReflection(pi.EVEN)
        system_one_combined.diagonalize()
        system_one_inverse = pi.SystemOne(system_one)
        system_one_inverse.setConservedParityUnderReflection(pi.ODD)
        system_one_inverse.diagonalize()
        system_one_combined.add(system_one_inverse)

        system_one.setConservedParityUnderReflection(pi.NA)
        system_one.diagonalize()

        # Build two atom system
        system_two = pi.SystemTwo(system_one, system_one, self.cache)
        system_two.setDistance(1)
        system_two.setOrder(3)
        system_two.diagonalize()

        # Note: it is important to use system_one_combined
        system_two_from_combined = pi.SystemTwo(system_one_combined, system_one_combined, self.cache)
        system_two_from_combined.setDistance(1)
        system_two_from_combined.setOrder(3)

        system_two_combined = None
        for sym_reflection, sym_inversion, sym_permutation in product([pi.EVEN, pi.ODD], [pi.EVEN, pi.ODD], [pi.EVEN, pi.ODD]):
            system_two_tmp = pi.SystemTwo(system_two_from_combined)
            system_two_tmp.setConservedParityUnderReflection(sym_reflection)
            system_two_tmp.setConservedParityUnderInversion(sym_inversion)
            system_two_tmp.setConservedParityUnderPermutation(sym_permutation)
            system_two_tmp.diagonalize()
            if system_two_combined is None: system_two_combined = system_two_tmp
            else: system_two_combined.add(system_two_tmp)
        system_two_combined.diagonalize()

        # Compare results
        w1 = np.sort(system_two_combined.diagonal)
        w2 = np.sort(system_two.diagonal)

        maxdiff = np.abs((w1-w2)/(np.max([w1,w2],axis=0)))
        maxdiff[np.abs(w1-w2)<1e-14] = np.abs(w1-w2)[np.abs(w1-w2)<1e-14]
        maxdiff = np.max(maxdiff)
        print("Two-atom system with combined binary symmetries, relative maximum deviation: ", maxdiff)
        self.assertAlmostEqual(maxdiff, 0, places=9)

    #######################################################
    ### Clean up ##########################################
    #######################################################

    def tearDown(self):
        # Delete cache object
        del self.cache
        
        # Delete cache directory
        shutil.rmtree(self.path_cache)

if __name__ == '__main__':
    unittest.main()
