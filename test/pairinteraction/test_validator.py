import json
import shutil
import unittest
from copy import deepcopy
from pathlib import Path

import pytest

from pairinteraction.validator.model import Model

directory = Path(__file__).parent
settings_directory = directory / "settings"
all_json_files = [f.stem for f in settings_directory.glob("*.json")]


class SimulationTests(unittest.TestCase):
    def testAllSettings(self):
        for name in all_json_files:
            self._oneTestModel(name)

    def _oneTestModel(self, name):
        self.loadSettings(name)
        output = self.settings.model_dump(exclude_unset=True)

        reference_path = directory / "data" / f"{name}__settings.json"
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
        self.removeCache()

    def testAtomErrors(self):
        with open(settings_directory / "atom1.json", encoding="utf-8") as f:
            original_settings = json.load(f)

        # Lists must have same length/steps
        settings = deepcopy(original_settings)
        settings["atom1"]["bfield_y"] = [0, 0.1, 0.2, 0.4]
        with pytest.raises(ValueError) as e_info:
            model = Model.model_validate(settings)
        assert "steps" in str(e_info.value)

        # States must have same species
        settings = deepcopy(original_settings)
        settings["atom1"]["states_of_interest"][0]["species"] = "Cs"
        with pytest.raises(ValueError) as e_info:
            model = Model.model_validate(settings)
        assert "species" in str(e_info.value)

        # Tests for is_real
        settings = deepcopy(original_settings)
        settings["atom1"]["efield_y"] = 0
        settings["atom1"]["bfield_y"] = [0, 0.1, 0.2]

        model = Model.model_validate(settings)
        assert not model.atom1.is_real

        settings["atom1"]["bfield_y"] = 0
        model = Model.model_validate(settings)
        assert model.atom1.is_real

    def loadSettings(self, name):
        with open(settings_directory / f"{name}.json", encoding="utf-8") as f:
            settings = json.load(f)
        self.settings = Model.model_validate(settings)

    def removeCache(self):
        if self.settings.numerics.path_cache is not None:
            shutil.rmtree(self.settings.numerics.path_cache, ignore_errors=True)
        self.settings = None


def update_settings_recursive(settings, new):
    for k, v in new.items():
        if isinstance(v, dict):
            settings[k] = update_settings_recursive(settings.get(k, {}), v)
        else:
            settings[k] = v
    return settings


if __name__ == "__main__":
    unittest.main()
