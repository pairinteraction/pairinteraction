import unittest
import sqlite3
import numpy as np
from multiprocessing import Pool
from pyinteraction import pireal as pi
#from calc import pairinteraction_real as pi

# TODO create and delete cache directory
# TODO calculate potential lines
# TODO check results

# TODO if the following is placed inside test_comparison, AttributeError: Can't pickle local object 'TestQuantumDefect.test_comparison.<locals>.calcEnergies'

# Define states
state_one = pi.StateOne("Rb", 61, 2, 1.5, 1.5)
state_two = pi.StateTwo(state_one, state_one)

#######################################################
### Check parallelization of one atom calculations ####
#######################################################

# Build one atom system
system_one = pi.SystemOne(state_one.element)
system_one.restrictEnergy(state_one.energy-40, state_one.energy+40)
system_one.restrictN(state_one.n-1, state_one.n+1)
system_one.restrictL(state_one.l-1, state_one.l+1)
system_one.setEfield([0, 0, 1])
system_one.buildInteraction() # will save time in the following

# Calculate stark shifts in parallel
efields = np.linspace(0,5,5)

def calcEnergies(field):
    tmp = pi.SystemOne(system_one) # copy allows to restrict energy further
    tmp.setEfield([0, 0, field])
    tmp.diagonalize()
    tmp.restrictEnergy(state_one.energy-20, state_one.energy+20)
    tmp.buildHamiltonian() # has to be called to apply the new restriction in energy
    return tmp

p = Pool()
stark_shifted_systems = np.array(p.map(calcEnergies, efields))
p.close()

#######################################################
### Check parallelization of two atom calculations ####
#######################################################

# Build two atom system
system_two = pi.SystemTwo(system_one, system_one)
system_two.restrictEnergy(state_two.energy-2, state_two.energy+2)
system_two.setConservedParityUnderPermutation(pi.ODD)
system_two.setDistance(1)
system_two.buildInteraction() # will save time in the following

# Calculate dipole interaction in parallel
distances = np.linspace(10,4,5)

def calcEnergies(distance):
    tmp = pi.SystemTwo(system_two) # copy allows to restrict energy further
    tmp.setDistance(distance)
    tmp.diagonalize()
    tmp.restrictEnergy(state_two.energy-0.1, state_two.energy+0.1)
    tmp.buildHamiltonian() # has to be called to apply the new restriction in energy
    return tmp

p = Pool()
dipole_interacting_systems = np.array(p.map(calcEnergies, distances))
p.close()

