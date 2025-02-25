#!/usr/bin/env -S uv run --script
# /// script
# requires-python = "==3.9"
# dependencies = [
# "pairinteraction==0.9.9",
# "numpy",
# ]
# ///

import json
import sys

import numpy as np
from timer import timer

from pairinteraction import pireal as pi


def main() -> None:
    settings = json.loads(sys.argv[1])
    species = settings["species"]
    n = settings["n"]
    l = settings["l"]
    j = settings["j"]
    m = settings["m"]
    delta_n = settings["delta_n"]
    delta_l = settings["delta_l"]
    delta_energy = settings["delta_energy_in_GHz"]
    distances = np.array(settings["distances_in_um"])
    order = settings["order"]

    with timer() as get_duration:
        cache = pi.MatrixElementCache()
        state_one = pi.StateOne(species, n, l, j, m)
        system_one = pi.SystemOne(state_one.getSpecies(), cache)
        system_one.restrictN(state_one.getN() - delta_n, state_one.getN() + delta_n)
        system_one.restrictL(state_one.getL() - delta_l, state_one.getL() + delta_l)
        state_two = pi.StateTwo(state_one, state_one)
        system_two = pi.SystemTwo(system_one, system_one, cache)
        system_two.restrictEnergy(state_two.getEnergy() - delta_energy, state_two.getEnergy() + delta_energy)
        system_two.setConservedMomentaUnderRotation([int(2 * m)])
        system_two.setConservedParityUnderInversion(pi.ODD)
        system_two.setConservedParityUnderPermutation(pi.ODD)
        system_two.setOrder(order)
        system_two.setDistance(np.min(distances))
        system_two.buildHamiltonian()
    duration_construction = get_duration()
    number_of_kets = system_two.getNumBasisvectors()
    number_of_states = number_of_kets

    with timer() as get_duration:
        for distance in distances:
            system_two.setDistance(distance)
            system_two.diagonalize()
    duration_diagonalization = get_duration()

    output = {
        "number_of_kets": number_of_kets,
        "number_of_states": number_of_states,
        "duration_construction": duration_construction,
        "duration_diagonalization": duration_diagonalization,
    }
    print("\n" + json.dumps(output))


if __name__ == "__main__":
    main()
