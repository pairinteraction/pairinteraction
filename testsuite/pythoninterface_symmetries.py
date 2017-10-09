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

        momenta = np.arange(-(state_one.j+1), (state_one.j+1)+1)

        # Diagonalize blockwise
        system_one_combined = None
        for m in momenta:
            if system_one_combined is None:
                system_one_combined = pi.SystemOne(system_one)
                system_one_combined.setConservedMomentaUnderRotation([m])
                system_one_combined.diagonalize()
            else:
                system_one_tmp = pi.SystemOne(system_one)
                system_one_tmp.setConservedMomentaUnderRotation([m])
                system_one_tmp.diagonalize()
                system_one_combined.add(system_one_tmp)

        # Diagonalize altogether
        system_one_alternative = pi.SystemOne(system_one)
        system_one_alternative.setConservedMomentaUnderRotation(momenta)
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
        system_two_even.diagonalize()

        system_two_odd = pi.SystemTwo(system_one_even, system_one_odd, path_cache)
        system_two_odd.add(pi.SystemTwo(system_one_odd, system_one_even, path_cache))
        system_two_odd.setDistance(1)
        system_two_odd.diagonalize()

        system_two_combined =  pi.SystemTwo(system_two_even)
        system_two_combined.add(system_two_odd)

        # Diagonalize altogether
        system_two = pi.SystemTwo(system_one, system_one, path_cache)
        system_two.setDistance(1)
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
        ### Clean up ##########################################
        #######################################################

        # Delete cache directory
        shutil.rmtree(path_cache)

if __name__ == '__main__':
    unittest.main()
