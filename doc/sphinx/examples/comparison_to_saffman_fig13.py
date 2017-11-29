import os
import numpy as np
import matplotlib.pyplot as plt
from multiprocessing import Pool
import sys
if sys.platform == "darwin": sys.path.append("/Applications/pairinteraction.app/Contents/Resources")
elif sys.platform == "win32": sys.path.append("C:\Program Files\pairinteraction")
from libpairinteraction import pireal as pi

distance = 10 # µm
bfields = np.linspace(0, 20, 200) # Gauss

state_one = pi.StateOne("Rb", 43, 2, 2.5, 0.5)
state_two = pi.StateTwo(state_one, state_one)

if not os.path.exists("./cache"):
    os.makedirs("./cache")

def setup_system_one(bfield):
    system_one = pi.SystemOne(state_one.element, "./cache")
    system_one.restrictEnergy(state_one.energy-100, state_one.energy+100)
    system_one.restrictN(state_one.n-2, state_one.n+2)
    system_one.restrictL(state_one.l-2, state_one.l+2)
    system_one.setBfield([0, 0, bfield])
    return system_one

def setup_system_two(system_one,angle):
    system_two = pi.SystemTwo(system_one, system_one, "./cache")
    system_two.restrictEnergy(state_two.energy-5, state_two.energy+5)
    system_two.setDistance(10)
    system_two.setAngle(angle)
    if angle == 0: system_two.setConservedMomentaUnderRotation([int(2*state_one.m)])
    system_two.setConservedParityUnderInversion(pi.ODD)
    system_two.setConservedParityUnderPermutation(pi.ODD)
    return system_two

def getEnergies(bfield):
    # Set up one atom system
    system_one = setup_system_one(bfield)
    system_one.diagonalize()

    # Calculate Zeeman shift
    zeemanshift = 2*system_one.diagonal[system_one.getVectorindex(state_one)]
    
    # Set up two atom system
    system_two = setup_system_two(system_one,angle)
    system_two.diagonalize()
    
    # Calculate blockade interaction
    eigenenergies = (system_two.diagonal-zeemanshift)*1e3 # MHz
    overlaps = system_two.getOverlap(state_two)
    blockade = 1/np.sqrt(np.sum(overlaps/eigenenergies**2))
    
    return blockade

if __name__ ==  '__main__':
    plt.xlabel(r"$B$ (Gauss)")
    plt.ylabel(r"Blockade (MHz)")
    plt.xlim(-1,20)
    plt.ylim(0,0.4)

    angle = 0 # rad
    with Pool() as pool:
        energies1 = pool.map(getEnergies, bfields)

    angle = np.pi/2 # rad
    with Pool() as pool:
        energies2 = pool.map(getEnergies, bfields)

    plt.plot(bfields, energies1, '-', label=r"$\theta = 0$")
    plt.plot(bfields, energies2, '-', label=r"$\theta = \pi/2$")
    plt.legend(loc=2, bbox_to_anchor=(1.02, 1), borderaxespad=0);