# /// script
# requires-python = "==3.13.*"
# dependencies = [
# "arc-alkali-rydberg-calculator",
# "numpy",
# ]
# ///

import json
import sys

import arc
import numpy as np
from timer import timer


def main() -> None:
    settings = json.loads(sys.argv[1])
    species = settings["species"]
    n = settings["n"]
    l = settings["l"]
    j = settings["j"]
    m = settings["m"]
    delta_n = settings["delta_n"]
    delta_l = settings["delta_l"]
    delta_energy = settings["delta_energy_in_GHz"] * 1e9
    distances = np.array(settings["distances_in_um"])
    order = settings["order"]

    assert species == "Rb", "Only Rb is supported."

    with timer() as get_duration:
        calc = arc.PairStateInteractions(arc.Rubidium(), n, l, j, n, l, j, m, m, interactionsUpTo=(order - 1) // 2)
        calc.defineBasis(0, 0, delta_n, delta_l, delta_energy, progressOutput=False)
    number_of_kets = len(calc.basisStates)
    duration_construction = get_duration()

    number_of_states = 50
    with timer() as get_duration:
        calc.diagonalise(distances, number_of_states, progressOutput=False)
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
