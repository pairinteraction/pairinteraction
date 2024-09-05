# https://pytest-qt.readthedocs.io/en/latest/tutorial.html
import io
import json
import unittest
import zipfile
from pathlib import Path

import numpy as np
from PyQt5.QtCore import Qt
from PyQt5.QtTest import QTest
from scipy.io import loadmat

from pairinteraction.gui.app import MainWindow

DATA_DIRECTORY = Path(__file__).parent / "data"
SKIP_CALC = True


def test_field_calc_button(qtbot, tmp_path):
    """Testing simulation single atom with all E and B fields on"""
    window = setup_window(qtbot, tmp_path)
    window.loadSettingsSystem(DATA_DIRECTORY / "field" / "settings.sconf")
    window.loadSettingsPlotter(DATA_DIRECTORY / "field" / "settings.pconf")
    calc_and_compare_energies(window, tmp_path, "field", dE=3)
    window.cleanupProcesses()


def test_potential_calc_button(qtbot, tmp_path):
    """Testing simulation for pairpotential"""
    window = setup_window(qtbot, tmp_path)
    window.loadSettingsSystem(DATA_DIRECTORY / "potential" / "settings.sconf")
    window.loadSettingsPlotter(DATA_DIRECTORY / "potential" / "settings.pconf")
    calc_and_compare_energies(window, tmp_path, "potential", dE=0.3)
    window.cleanupProcesses()


def setup_window(qtbot, tmp_path):
    window = MainWindow()
    window.show()
    qtbot.addWidget(window)
    window.path_cache = tmp_path / "cache"
    window.ui.action_sconf_reset.trigger()
    window.ui.action_pconf_reset.trigger()
    return window


def calc_and_compare_energies(window, tmp_path, ref_data, dE, dE_tol=1e-3, use_python_api="default"):
    if use_python_api == "default":
        for use_python_api in [False, True]:
            calc_and_compare_energies(window, tmp_path, ref_data, dE, use_python_api=use_python_api)
        return
    window.ui.checkbox_use_python_api.setChecked(use_python_api)
    window.autosetSymmetrization()

    if ref_data == "field":
        widget_calc = window.ui.pushbutton_field1_calc
        widget_save = window.ui.pushbutton_field1_save
    elif ref_data == "potential":
        widget_calc = window.ui.pushbutton_potential_calc
        widget_save = window.ui.pushbutton_potential_save

    # Run simulation
    if not SKIP_CALC:
        QTest.mouseClick(widget_calc, Qt.LeftButton)
    while window.timer.isActive():
        window.checkForData()

    # Save current data
    window.savePlot = False
    tmp_data_path = tmp_path / "data"
    tmp_data_path.mkdir(exist_ok=True)
    forceFilename = tmp_data_path / "current"
    window.forceFilename = forceFilename
    QTest.mouseClick(widget_save, Qt.LeftButton)

    data = {}
    sconfig = {}
    # Load current data
    with zipfile.ZipFile(forceFilename, "r") as zip_file:
        with zip_file.open("data.mat") as f:
            f_io = io.BytesIO(f.read())
            data["current"] = loadmat(f_io)
        with zip_file.open("settings.sconf") as f:
            sconfig["current"] = json.load(f)

    # Load reference data
    with open(DATA_DIRECTORY / ref_data / "data.mat", "rb") as f:
        data["ref"] = loadmat(f)
    with open(DATA_DIRECTORY / ref_data / "settings.sconf") as f:
        sconfig["ref"] = json.load(f)

    if SKIP_CALC:
        return

    # Check if configs match # unecessary since we load the same config
    for k, v in sconfig["ref"].items():
        assert sconfig["current"][k] == v

    if len(data["current"]["eigenvalues"]) == 0:
        raise ValueError("No eigenvalues found in current data, some bug in calculating the eigenvalues.")

    # Check if central eigenvalues (+/- dE) match
    for i in range(len(data["ref"]["eigenvalues"])):
        Es = {k: np.array(mat["eigenvalues"])[i] for k, mat in data.items()}
        Es = {k: E[np.abs(E) < dE] for k, E in Es.items()}
        assert len(np.abs(Es["ref"]) == len(Es["current"]))
        diff_rel = np.abs(Es["ref"] - Es["current"])
        assert np.all(diff_rel <= dE_tol)


if __name__ == "__main__":
    unittest.main()
