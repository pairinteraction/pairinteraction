import shutil
import unittest
from pathlib import Path

import numpy as np

from pairinteraction.model.simulation import ModelSimulation
from pairinteraction.simulation_old.simulation import Simulation

directory = Path(__file__).parent
names = ["simulation_1", "simulation_2", "simulation_3"]


class SimulationTests(unittest.TestCase):
    def test_simulations(self):
        for name in names:
            self.one_test_simulation(name)

    def one_test_simulation(self, name):
        self.model = ModelSimulation.model_validate_json_file(directory / "models" / f"{name}.json")
        simulation = Simulation(self.model)
        results_list = simulation.run()
        energies_list = [results["energies"] for results in results_list]

        reference_path = directory / "data" / f"{name}__energies.txt"
        # use this for updating the reference
        # self.save_energies_list(reference_path, energies_list)

        reference_list = self.load_energies_list(reference_path)

        assert len(energies_list) == len(reference_list)
        for energies, reference in zip(energies_list, reference_list):
            assert len(energies) == len(reference)
            assert np.allclose(energies, reference, rtol=1e-4)

        self.remove_cache()

    def save_energies_list(self, name, energies_list):
        with open(name, "w", encoding="utf-8") as f:
            for sublist in energies_list:
                f.write(" ".join(map(str, sublist)) + "\n")

    def load_energies_list(self, name):
        with open(name, encoding="utf-8") as f:
            return [list(map(float, line.split())) for line in f]

    def remove_cache(self):
        if self.model.numerics.path_cache is not None:
            shutil.rmtree(self.model.numerics.path_cache, ignore_errors=True)
        self.model = None


if __name__ == "__main__":
    unittest.main()
