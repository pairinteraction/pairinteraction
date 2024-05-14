import json
import shutil
import unittest
from copy import deepcopy
from pathlib import Path

from pairinteraction.model.model import ModelSimulation

directory = Path(__file__).parent
models_directory = directory / "models"
all_json_files = [f.stem for f in models_directory.glob("*.json")]


class SimulationTests(unittest.TestCase):
    def test_all_models(self):
        for name in all_json_files:
            self._one_test_model(name)

    def _one_test_model(self, name):
        self.model = ModelSimulation.model_validate_json_file(models_directory / f"{name}.json")
        output = self.model.model_dump(exclude_unset=True)

        reference_path = directory / "data" / f"{name}__model.json"
        # use this for updating the reference
        # with open(reference_path, "w", encoding="utf-8") as f:
        #     json.dump(output, f, indent=4)
        with open(reference_path, encoding="utf-8") as f:
            reference = json.load(f)

        def assert_dict_equal(d1, d2):
            assert len(d1) == len(d2)
            for k, v in d1.items():
                if isinstance(v, dict):
                    assert_dict_equal(v, d2[k])
                else:
                    assert d2[k] == v, f"{d2[k]} != {v}"

        assert_dict_equal(reference, output)
        self.remove_cache()

    def test_atom_errors(self):
        with open(models_directory / "atom1.json", encoding="utf-8") as f:
            original_model_dict = json.load(f)

        # Lists must have same length/steps
        model = deepcopy(original_model_dict)
        model["atom1"]["bfield_y"] = [0, 0.1, 0.2, 0.4]
        with self.assertRaises(ValueError) as cm:
            model = ModelSimulation.model_validate(model)
        assert "steps" in str(cm.exception)

        # States must have same species
        model = deepcopy(original_model_dict)
        model["atom1"]["states_of_interest"][0]["species"] = "Cs"
        with self.assertRaises(ValueError) as cm:
            model = ModelSimulation.model_validate(model)
        assert "species" in str(cm.exception)

        # Tests for is_real
        model = deepcopy(original_model_dict)
        model["atom1"]["efield_y"] = 0
        model["atom1"]["bfield_y"] = [0, 0.1, 0.2]
        model = ModelSimulation.model_validate(model)
        assert not model.atom1.is_real

        model = deepcopy(original_model_dict)
        model["atom1"]["efield_y"] = 0
        model["atom1"]["bfield_y"] = 0
        model = ModelSimulation.model_validate(model)
        assert model.atom1.is_real

    def remove_cache(self):
        if self.model.numerics.path_cache is not None:
            shutil.rmtree(self.model.numerics.path_cache, ignore_errors=True)
        self.model = None


if __name__ == "__main__":
    unittest.main()
