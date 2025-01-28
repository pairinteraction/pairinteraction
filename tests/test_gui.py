# https://pytest-qt.readthedocs.io/en/latest/tutorial.html
import io
import json
import zipfile
from pathlib import Path

import numpy as np
import pytest
from PyQt5.QtCore import Qt
from PyQt5.QtTest import QTest
from scipy.io import loadmat

from pairinteraction.gui.app import MainWindow

reference_directory = Path(__file__).parent.parent / "data/reference_gui"

pytest.skip("Skip gui tests for now", allow_module_level=True)


def setup_window(qtbot, tmp_path) -> MainWindow:
    window = MainWindow()
    qtbot.addWidget(window)
    window.path_cache = tmp_path / "cache"
    window.ui.action_sconf_reset.trigger()
    window.ui.action_pconf_reset.trigger()
    return window


def calc_and_compare_energies(generate_reference, window, tmp_path, ref_data, use_python_api, dE, dE_tol=1e-3) -> None:
    window.ui.checkbox_use_python_api.setChecked(use_python_api)
    window.autosetSymmetrization()

    if ref_data == "field":
        widget_calc = window.ui.pushbutton_field1_calc
        widget_save = window.ui.pushbutton_field1_save
    elif ref_data == "potential":
        widget_calc = window.ui.pushbutton_potential_calc
        widget_save = window.ui.pushbutton_potential_save

    reference_data_file = reference_directory / ref_data / "data.mat"
    reference_settings_file = reference_directory / ref_data / "settings.sconf"

    pytest.skip("Calculations using the GUI are not yet supported.")

    # Run simulation
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

    if generate_reference:
        with zipfile.ZipFile(forceFilename, "r") as zip_file:
            with zip_file.open("data.mat") as f:
                with open(reference_data_file, "wb") as f_ref:
                    f_ref.write(f.read())
            with zip_file.open("settings.sconf") as f:
                with open(reference_settings_file, "w") as f_ref:
                    f_ref.write(f.read())
        pytest.skip("Reference data generated, skipping comparison test")

    # Load current and reference data
    data = {}
    sconfig = {}

    with zipfile.ZipFile(forceFilename, "r") as zip_file:
        with zip_file.open("data.mat") as f:
            f_io = io.BytesIO(f.read())
            data["current"] = loadmat(f_io)
        with zip_file.open("settings.sconf") as f:
            sconfig["current"] = json.load(f)

    with open(reference_data_file, "rb") as f:
        data["ref"] = loadmat(f)
    with open(reference_settings_file) as f:
        sconfig["ref"] = json.load(f)

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


@pytest.mark.parametrize("use_python_api", [False, True])
def test_field_calc_button(generate_reference, qtbot, tmp_path, use_python_api) -> None:
    """Testing simulation single atom with all E and B fields on."""
    window = setup_window(qtbot, tmp_path)
    window.loadSettingsSystem(reference_directory / "field" / "settings.sconf")
    window.loadSettingsPlotter(reference_directory / "field" / "settings.pconf")
    calc_and_compare_energies(generate_reference, window, tmp_path, "field", use_python_api, dE=3)
    window.cleanupProcesses()


@pytest.mark.parametrize("use_python_api", [False, True])
def test_potential_calc_button(generate_reference, qtbot, tmp_path, use_python_api) -> None:
    """Testing simulation for pairpotential."""
    window = setup_window(qtbot, tmp_path)
    window.loadSettingsSystem(reference_directory / "potential" / "settings.sconf")
    window.loadSettingsPlotter(reference_directory / "potential" / "settings.pconf")
    calc_and_compare_energies(generate_reference, window, tmp_path, "potential", use_python_api, dE=0.3)
    window.cleanupProcesses()
