import unittest
import sqlite3
import numpy as np
from scipy.sparse import coo_matrix
from multiprocessing import Pool
from @LIBNAME@ import pireal as pi
import tempfile
import shutil

class TestPythoninterfaceMultiprocessing(unittest.TestCase):

    def test_comparison(self):
        
        #######################################################
        ### Preparations ######################################
        #######################################################

        # Create cache directory
        path_cache = tempfile.mkdtemp()
        
        #######################################################
        ### Check parallelization of one atom calculations ####
        #######################################################

        # Define state
        state_one = pi.StateOne("Rb", 61, 2, 1.5, 1.5)
        
        # Build one atom system
        system_one = pi.SystemOne(state_one.element, path_cache)
        system_one.restrictEnergy(state_one.energy-40, state_one.energy+40)
        system_one.restrictN(state_one.n-1, state_one.n+1)
        system_one.restrictL(state_one.l-1, state_one.l+1)
        system_one.setEfield([0, 0, 1])
        system_one.buildInteraction() # will save time in the following

        # Calculate stark shifts in parallel
        efields = np.linspace(0,5,5)

        global calcEnergies_one
        def calcEnergies_one(field):
            tmp = pi.SystemOne(system_one) # copy allows to restrict energy further
            tmp.setEfield([0, 0, field])
            tmp.diagonalize()
            tmp.restrictEnergy(state_one.energy-20, state_one.energy+20)
            tmp.buildHamiltonian() # has to be called to apply the new restriction in energy
            return tmp

        p = Pool(2)
        stark_shifted_systems = np.array(p.map(calcEnergies_one, efields))
        p.close()

        # Get potential lines
        line_mapper = np.array([])
        line_idx = []
        line_step = []
        line_val = []

        for i in range(len(stark_shifted_systems)):

            # Get line segments pointing from iFirst to iSecond
            if i > 0:
                connections = stark_shifted_systems[i-1].getConnections(stark_shifted_systems[i], 0.001)
                iFirst = np.array(connections[0])
                iSecond = np.array(connections[1])

            else:
                iFirst = np.arange(stark_shifted_systems[i].getNumVectors())
                iSecond = iFirst

            if len(iFirst) > 0:
                # Enlarge the mapper of line indices so that every iFirst is mapped to a line index
                line_mapper.resize(np.max(iFirst)+1)
                boolarr = line_mapper[iFirst] == 0
                tmp = line_mapper[iFirst]
                tmp[boolarr] = np.max(line_mapper)+1+np.arange(np.sum(boolarr))
                line_mapper[iFirst] = tmp

                # Store line segments    
                line_idx += line_mapper[iFirst].tolist()
                line_step += (i*np.ones(len(iFirst))).tolist()
                line_val += stark_shifted_systems[i].diagonal[iSecond].tolist()

                # Assign the line indices to iSecond, delete line indices that could not be assigned
                tmp = line_mapper[iFirst]
                line_mapper = np.zeros(np.max(iSecond)+1)
                line_mapper[iSecond] = tmp
        
        lines_energies = coo_matrix((line_val, (np.array(line_idx)-1,line_step))).toarray()
        lines_energies[lines_energies == 0] = np.nan

        self.assertEqual(len(lines_energies), 30)
        self.assertEqual(np.sum(np.isnan(lines_energies)), 2)

        #######################################################
        ### Check parallelization of two atom calculations ####
        #######################################################

        # Define state
        state_two = pi.StateTwo(state_one, state_one)

        # Build two atom system
        system_two = pi.SystemTwo(system_one, system_one, path_cache)
        system_two.restrictEnergy(state_two.energy-2, state_two.energy+2)
        system_two.setConservedParityUnderPermutation(pi.ODD)
        system_two.setDistance(1)
        system_two.buildInteraction() # will save time in the following

        # Calculate dipole interaction in parallel
        distances = np.linspace(10,4,5)

        global calcEnergies_two
        def calcEnergies_two(distance):
            tmp = pi.SystemTwo(system_two) # copy allows to restrict energy further
            tmp.setDistance(distance)
            tmp.diagonalize()
            tmp.restrictEnergy(state_two.energy-0.1, state_two.energy+0.1)
            tmp.buildHamiltonian() # has to be called to apply the new restriction in energy
            return tmp

        p = Pool(2)
        dipole_interacting_systems = np.array(p.map(calcEnergies_two, distances))
        p.close()

        # Get potential lines
        line_mapper = np.array([])
        line_idx = []
        line_step = []
        line_val = []

        for i in range(len(dipole_interacting_systems)):

            # Get line segments pointing from iFirst to iSecond
            if i > 0:
                connections = dipole_interacting_systems[i-1].getConnections(dipole_interacting_systems[i], 0.001)
                iFirst = np.array(connections[0])
                iSecond = np.array(connections[1])

            else:
                iFirst = np.arange(dipole_interacting_systems[i].getNumVectors())
                iSecond = iFirst

            if len(iFirst) > 0:
                # Enlarge the mapper of line indices so that every iFirst is mapped to a line index
                line_mapper.resize(np.max(iFirst)+1)
                boolarr = line_mapper[iFirst] == 0
                tmp = line_mapper[iFirst]
                tmp[boolarr] = np.max(line_mapper)+1+np.arange(np.sum(boolarr))
                line_mapper[iFirst] = tmp

                # Store line segments    
                line_idx += line_mapper[iFirst].tolist()
                line_step += (i*np.ones(len(iFirst))).tolist()
                line_val += dipole_interacting_systems[i].diagonal[iSecond].tolist()

                # Assign the line indices to iSecond, delete line indices that could not be assigned
                tmp = line_mapper[iFirst]
                line_mapper = np.zeros(np.max(iSecond)+1)
                line_mapper[iSecond] = tmp

        lines_energies = coo_matrix((line_val, (np.array(line_idx)-1,line_step))).toarray()
        lines_energies[lines_energies == 0] = np.nan

        self.assertEqual(len(lines_energies), 10)
        self.assertEqual(np.sum(np.isnan(lines_energies)), 0)
        
        #######################################################
        ### Clean up ##########################################
        #######################################################

        # Delete cache directory
        shutil.rmtree(path_cache)

if __name__ == '__main__':
    unittest.main()
