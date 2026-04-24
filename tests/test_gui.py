# SPDX-FileCopyrightText: 2024 PairInteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

from __future__ import annotations

from typing import TYPE_CHECKING, Any, Literal

import numpy as np
import pytest
from pairinteraction_gui.calculate.calculate_one_atom import _get_surface_position
from pairinteraction_gui.main_window import MainWindow
from pairinteraction_gui.page.one_atom_page import OneAtomPage

from .utils import REFERENCE_PATHS, compare_eigensystem_to_reference

if TYPE_CHECKING:
    from pathlib import Path

    from pairinteraction_gui.page import LifetimesPage, OneAtomPage
    from pairinteraction_gui.page.two_atoms_page import TwoAtomsPage
    from pytestqt.qtbot import QtBot


@pytest.fixture
def base_window(qtbot: QtBot, tmp_path: Path) -> MainWindow:
    window = MainWindow(cache_dir=tmp_path)
    window.show()
    qtbot.addWidget(window)
    return window


@pytest.fixture
def window_starkmap(base_window: MainWindow) -> MainWindow:
    one_atom_page: OneAtomPage = base_window.stacked_pages.getNamedWidget("OneAtomPage")  # type: ignore [assignment]
    one_atom_page.ket_config.species_combo_list[0].setCurrentText("Rb")
    ket_qn = one_atom_page.ket_config.stacked_qn_list[0].currentWidget()
    ket_qn.items["n"].setValue(60)
    ket_qn.items["l"].setValue(0)
    ket_qn.items["m"].setValue(0.5)

    basis_qn = one_atom_page.basis_config.stacked_basis_list[0].currentWidget()
    basis_qn.items["n"].setValue(2)
    basis_qn.items["l"].setValue(2)
    basis_qn.items["m"].setChecked(False)

    calculation_config = one_atom_page.calculation_config
    calculation_config.steps.setValue(11)
    system_config = one_atom_page.system_config
    system_config.Ez.spinboxes[1].setValue(10)

    return base_window


@pytest.fixture
def window_pair_potential(base_window: MainWindow) -> MainWindow:
    two_atoms_page: TwoAtomsPage = base_window.stacked_pages.getNamedWidget("TwoAtomsPage")  # type: ignore [assignment]
    two_atoms_page.ket_config.species_combo_list[0].setCurrentText("Rb")
    for ket_qn_stacked in two_atoms_page.ket_config.stacked_qn_list:
        ket_qn = ket_qn_stacked.currentWidget()
        ket_qn.items["n"].setValue(60)
        ket_qn.items["l"].setValue(0)
        ket_qn.items["m"].setValue(0.5)

    for basis_qn_stacked in two_atoms_page.basis_config.stacked_basis_list:
        basis_qn = basis_qn_stacked.currentWidget()
        basis_qn.items["n"].setValue(2)
        basis_qn.items["l"].setValue(2)
        basis_qn.items["m"].setChecked(False)

    two_atoms_page.basis_config.pair_delta_energy.setValue(3)
    two_atoms_page.basis_config.pair_m_range.setValues(1, 1)

    calculation_config = two_atoms_page.calculation_config
    calculation_config.steps.setValue(5)
    system_config = two_atoms_page.system_config
    system_config.distance.setValues(1, 5)
    return base_window


@pytest.fixture
def window_lifetimes(base_window: MainWindow) -> MainWindow:
    lifetimes_page: LifetimesPage = base_window.stacked_pages.getNamedWidget("LifetimesPage")  # type: ignore [assignment]
    lifetimes_page.ket_config.species_combo_list[0].setCurrentText("Rb")
    ket_qn = lifetimes_page.ket_config.stacked_qn_list[0].currentWidget()
    ket_qn.items["n"].setValue(60)
    ket_qn.items["l"].setValue(0)
    ket_qn.items["m"].setValue(0.5)
    lifetimes_page.ket_config.item_temperature.setValue(300)
    return base_window


def test_main_window_basic(qtbot: QtBot, window_starkmap: MainWindow) -> None:
    """Test basic main window functionality."""
    one_atom_page: OneAtomPage = window_starkmap.stacked_pages.getNamedWidget("OneAtomPage")  # type: ignore [assignment]
    qn_item = one_atom_page.ket_config.stacked_qn_list[0].currentWidget().items["n"]
    qn_item.setValue(60)

    ket_label = one_atom_page.ket_config.ket_label_list[0].text()
    assert all(x in ket_label for x in ["Rb", "60", "S", "1/2"])
    assert qn_item.label.text() == "n"
    assert qn_item.value() == 60

    qn_item.setValue(61)
    ket_label = one_atom_page.ket_config.ket_label_list[0].text()
    assert qn_item.value() == 61
    assert all(x in ket_label for x in ["Rb", "61", "S", "1/2"])

    # make the basis smaller for faster test
    basis_qn = one_atom_page.basis_config.stacked_basis_list[0].currentWidget()
    basis_qn.items["n"].setValue(1)
    basis_qn.items["l"].setValue(1)
    basis_qn.items["m"].setValue(0)

    one_atom_page.calculate_and_abort.getNamedWidget("Calculate").click()
    qtbot.waitUntil(lambda: one_atom_page._calculation_finished, timeout=30_000)  # ci macOS-15-intel is very slow
    qtbot.waitUntil(lambda: one_atom_page._plot_finished, timeout=5_000)
    window_starkmap.close()


def test_one_atom_page(window_starkmap: MainWindow) -> None:
    _test_calculate_page(window_starkmap, "OneAtomPage", "stark_map")


def test_one_atom_surface_angle_range(base_window: MainWindow) -> None:
    page: OneAtomPage = base_window.stacked_pages.getNamedWidget("OneAtomPage")  # type: ignore [assignment]
    page.ket_config.species_combo_list[0].setCurrentText("Rb")
    ket_qn = page.ket_config.stacked_qn_list[0].currentWidget()
    ket_qn.items["n"].setValue(60)
    ket_qn.items["l"].setValue(0)
    ket_qn.items["m"].setValue(0.5)

    basis_qn = page.basis_config.stacked_basis_list[0].currentWidget()
    basis_qn.items["n"].setValue(2)
    basis_qn.items["l"].setValue(2)
    basis_qn.items["m"].setChecked(False)

    system_config = page.system_config
    calculation_config = page.calculation_config

    calculation_config.steps.setValue(3)
    assert system_config.distance_to_surface.label.text() == "Distance to surface"
    system_config.distance_to_surface.setChecked(True)
    system_config.distance_to_surface.setValues(5, 5)
    system_config.surface_angle.setChecked(True)
    system_config.surface_angle.setValues(0, 90)

    ranges = system_config.get_ranges_dict()
    assert ranges["DistanceToSurface"] == [5.0, 5.0, 5.0]
    assert ranges["SurfaceAngle"] == [0.0, 45.0, 90.0]

    python_code = page._create_python_code()
    assert "surface_angle = np.linspace(0.0, 90.0, steps)" in python_code
    assert "distance_to_surface = np.linspace(5.0, 5.0, steps)" in python_code
    assert "surface_position = distance_to_surface[i] * np.array(" in python_code
    assert "normal=[np.sin(surface_angle_rad), 0, np.cos(surface_angle_rad)]" in python_code


def test_one_atom_surface_distance_follows_surface_normal() -> None:
    assert np.allclose(_get_surface_position(5.0, 0.0), [0.0, 0.0, 5.0])
    assert np.allclose(_get_surface_position(5.0, 90.0), [5.0, 0.0, 0.0])
    assert np.allclose(_get_surface_position(5.0, 45.0), [5 / np.sqrt(2), 0.0, 5 / np.sqrt(2)])


def test_two_atoms_page(window_pair_potential: MainWindow) -> None:
    _test_calculate_page(window_pair_potential, "TwoAtomsPage", "pair_potential")


def test_lifetimes_page(window_lifetimes: MainWindow) -> None:
    page: LifetimesPage = window_lifetimes.stacked_pages.getNamedWidget("LifetimesPage")  # type: ignore [assignment]
    parameters, results = page.calculate()
    page.plotwidget.plot(parameters, results)

    assert results.lifetime > 0
    assert len(results.kets_sp) == len(results.transition_rates_sp)
    assert len(results.kets_bbr) == len(results.transition_rates_bbr)

    python_code = page._create_python_code()
    python_code = python_code.replace("plt.show()", "plt.close()")  # HACK, otherwise it will block the test

    locals_globals: dict[str, Any] = {}
    exec(python_code, locals_globals, locals_globals)  # noqa: S102

    rates_dict = {
        key: [sum(r for _, r in page.plotwidget.sorted_rates[label][n]) for n in locals_globals["n_list"]]
        for key, label in [("SP", "Spontaneous Decay"), ("BBR", "Black Body Radiation")]
    }

    assert np.isclose(locals_globals["lifetime"], results.lifetime)
    for key, rates in rates_dict.items():
        assert np.allclose(locals_globals["rates_summed"][key], rates)


def _test_calculate_page(
    window: MainWindow,
    page_name: Literal["OneAtomPage", "TwoAtomsPage"],
    reference_name: str,
) -> None:
    page: OneAtomPage | TwoAtomsPage = window.stacked_pages.getNamedWidget(page_name)  # type: ignore [assignment]
    ket_energy_0 = sum(page.ket_config.get_ket_atom(i).get_energy("GHz") for i in range(page.ket_config.n_atoms))

    # Test calculation with fast mode off
    page.calculation_config.fast_mode.setChecked(False)
    _parameters, results = page.calculate()
    energies = np.array(results.energies) + ket_energy_0
    compare_eigensystem_to_reference(REFERENCE_PATHS[reference_name], energies, np.array(results.ket_overlaps))

    # Test calculation with fast mode on and diagonalze energy_range
    # NOTE: with fast mode, the overlaps are different, so we don't compare them
    page.calculation_config.fast_mode.setChecked(True)
    page.calculation_config.energy_range.setChecked(True)
    page.calculation_config.energy_range.setValues(-200, 200)
    parameters, results = page.calculate()
    energies = np.array(results.energies) + ket_energy_0
    compare_eigensystem_to_reference(REFERENCE_PATHS[reference_name], energies)
    page.plotwidget.plot(parameters, results)

    if page_name == "TwoAtomsPage":
        # Test fitting of c3/c6/c3+c6 potential curves
        # We need to pass the data to the plotwidget.
        page.plotwidget.fit("c3")
        # We run a test fit twice, as it should iterate through the potential curves
        page.plotwidget.fit("c3")
        page.plotwidget.fit("c6")
        page.plotwidget.fit("c3+c6")
        # One could consider adding some return values to check if the estimated values are reasonable.
        # Currently fit does not return anything, it just plots and shows the fit parameters in UI

    # Test export to Python code
    python_code = page._create_python_code()
    python_code = python_code.replace("plt.show()", "plt.close()")  # HACK, otherwise it will block the test

    # HACK, see also https://stackoverflow.com/questions/45132645/list-comprehension-in-exec-with-empty-locals-nameerror
    locals_globals: dict[str, Any] = {}
    exec(python_code, locals_globals, locals_globals)  # noqa: S102
    energies = np.array(locals_globals["energies_list"]) + ket_energy_0
    compare_eigensystem_to_reference(REFERENCE_PATHS[reference_name], energies)


def test_save_and_restore_settings(qtbot: QtBot, tmp_path: Path) -> None:
    window = MainWindow(cache_dir=tmp_path)
    window.show()
    qtbot.addWidget(window)
    one_atom_page: OneAtomPage = window.stacked_pages.getNamedWidget("OneAtomPage")  # type: ignore [assignment]
    ket_qn = one_atom_page.ket_config.stacked_qn_list[0].currentWidget()
    ket_qn.items["n"].setValue(222)
    ket_qn.items["l"].setValue(999)
    window.close()

    window = MainWindow(cache_dir=tmp_path)
    window.show()
    qtbot.addWidget(window)
    one_atom_page = window.stacked_pages.getNamedWidget("OneAtomPage")  # type: ignore [assignment]
    ket_qn = one_atom_page.ket_config.stacked_qn_list[0].currentWidget()

    assert ket_qn.items["n"].value() == 222
    assert ket_qn.items["l"].value() == 999
