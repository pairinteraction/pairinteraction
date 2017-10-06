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
        ### Check reflection symmetry of one atom systems #####
        #######################################################

        # Define state
        state_one = pi.StateOne("Rb", 61, 1, 0.5, 0.5)

        # Build one atom system
        system_one = pi.SystemOne(state_one.element, path_cache)
        system_one.restrictEnergy(state_one.energy-40, state_one.energy+40)
        system_one.restrictN(state_one.n-1, state_one.n+1)
        system_one.restrictL(state_one.l-1, state_one.l+1)
        system_one.restrictJ(state_one.j-1, state_one.j+1)
        system_one.setBfield([0, 100, 0])
        system_one.setEfield([1, 0, 2])
                
        # Diagonalize blockwise
        system_one_combined = pi.SystemOne(system_one)
        system_one_combined.setConservedParityUnderReflection(pi.EVEN)
        system_one_combined.diagonalize()
                
        system_one_odd = pi.SystemOne(system_one)
        system_one_odd.setConservedParityUnderReflection(pi.ODD)
        system_one_odd.diagonalize()
        
        system_one_combined.add(system_one_odd)
                
        # Diagonalize altogether
        system_one.setConservedParityUnderReflection(pi.NA)
        system_one.diagonalize()
        
        # Compare results
        w1 = np.sort(system_one_combined.diagonal)
        w2 = np.sort(system_one.diagonal)
        
        maxdiff = np.max(np.abs((w1-w2)/(1e-64+np.min([w1,w2],axis=0)))) # TODO check without the need of +1e-64
        print("One-atom system with reflection symmetry, relative maximum deviation: ", maxdiff)
        self.assertAlmostEqual(maxdiff, 0, places=6)
                        
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
        
        maxdiff12 = np.max(np.abs((w1-w2)/(1e-64+np.min([w1,w2],axis=0)))) # TODO check without the need of +1e-64
        maxdiff23 = np.max(np.abs((w3-w2)/(1e-64+np.min([w3,w2],axis=0)))) # TODO check without the need of +1e-64
        print("One-atom system with rotation symmetry, relative maximum deviation: ", maxdiff12, ", ", maxdiff23)
        self.assertAlmostEqual(maxdiff12, 0, places=6)
        self.assertAlmostEqual(maxdiff23, 0, places=6)

        #######################################################
        ### Clean up ##########################################
        #######################################################

        # Delete cache directory
        shutil.rmtree(path_cache)

if __name__ == '__main__':
    unittest.main()
