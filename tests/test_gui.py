# SPDX-FileCopyrightText: 2024 Pairinteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

from typing import TYPE_CHECKING

import numpy as np
import pytest
from pairinteraction_gui.calculate.calculate_two_atoms import ParametersTwoAtoms, _calculate_two_atoms
from pairinteraction_gui.main_window import MainWindow

from .test_pair_potential import compare_pair_potential_to_reference
from .test_starkmap import compare_starkmap_to_reference

if TYPE_CHECKING:
    from pairinteraction_gui.page import OneAtomPage
    from pairinteraction_gui.page.two_atoms_page import TwoAtomsPage
    from pytestqt.qtbot import QtBot


@pytest.fixture
def window_starkmap(qtbot: "QtBot") -> MainWindow:
    window = MainWindow()
    window.show()
    qtbot.addWidget(window)

    one_atom_page: OneAtomPage = window.stacked_pages.getNamedWidget("OneAtomPage")
    one_atom_page.ket_config.species_combo[0].setCurrentText("Rb")
    ket_qn = one_atom_page.ket_config.stacked_qn[0].currentWidget()
    ket_qn.items["n"].setValue(60)
    ket_qn.items["l"].setValue(0)
    ket_qn.items["m"].setValue(0.5)

    basis_qn = one_atom_page.basis_config.stacked_basis[0].currentWidget()
    basis_qn.items["n"].setValue(2)
    basis_qn.items["l"].setValue(2)
    basis_qn.items["m"].setChecked(False)

    calculation_config = one_atom_page.calculation_config
    calculation_config.steps.setValue(11)
    system_config = one_atom_page.system_config
    system_config.Ez.spinboxes[1].setValue(10)

    return window


def test_main_window_basic(qtbot: "QtBot", window_starkmap: "MainWindow") -> None:
    """Test basic main window functionality."""
    one_atom_page: OneAtomPage = window_starkmap.stacked_pages.getNamedWidget("OneAtomPage")  # type: ignore [assignment]
    qn_item = one_atom_page.ket_config.stacked_qn[0].currentWidget().items["n"]
    qn_item.setValue(60)

    ket_label = one_atom_page.ket_config.ket_label[0].text()
    assert all(x in ket_label for x in ["Rb", "60", "S", "1/2"])
    assert qn_item.label.text() == "n"
    assert qn_item.value() == 60

    qn_item.setValue(61)
    ket_label = one_atom_page.ket_config.ket_label[0].text()
    assert qn_item.value() == 61
    assert all(x in ket_label for x in ["Rb", "61", "S", "1/2"])

    # make the basis smaller for faster test
    basis_qn = one_atom_page.basis_config.stacked_basis[0].currentWidget()
    basis_qn.items["n"].setValue(1)
    basis_qn.items["l"].setValue(1)
    basis_qn.items["m"].setValue(0)

    one_atom_page.calculate_and_abort.getNamedWidget("Calculate").click()
    qtbot.waitUntil(lambda: one_atom_page._calculation_finished, timeout=20_000)  # ci macOS-13 is very slow
    qtbot.waitUntil(lambda: one_atom_page._plot_finished, timeout=5_000)
    window_starkmap.close()


def test_calculate_one_atom(window_starkmap: "MainWindow") -> None:
    one_atom_page: OneAtomPage = window_starkmap.stacked_pages.getNamedWidget("OneAtomPage")  # type: ignore [assignment]
    ket = one_atom_page.ket_config.get_ket_atom(0)

    one_atom_page.calculation_config.fast_mode.setChecked(False)
    _parameters, results = one_atom_page.calculate()
    energies = np.array(results.energies) + ket.get_energy("GHz")
    compare_starkmap_to_reference(energies, np.array(results.ket_overlaps))

    one_atom_page.calculation_config.fast_mode.setChecked(True)
    _parameters, results = one_atom_page.calculate()
    energies = np.array(results.energies) + ket.get_energy("GHz")
    compare_starkmap_to_reference(energies)  # with fast mode, the overlaps are different, so we don't compare them


def test_calculate_two_atoms(qtbot: "QtBot") -> None:
    window = MainWindow()
    window.show()
    qtbot.addWidget(window)

    two_atoms_page: TwoAtomsPage = window.stacked_pages.getNamedWidget("TwoAtomsPage")
    two_atoms_page.ket_config.species_combo[0].setCurrentText("Rb")
    for ket_qn_stacked in two_atoms_page.ket_config.stacked_qn:
        ket_qn = ket_qn_stacked.currentWidget()
        ket_qn.items["n"].setValue(60)
        ket_qn.items["l"].setValue(0)
        ket_qn.items["m"].setValue(0.5)

    for basis_qn_stacked in two_atoms_page.basis_config.stacked_basis:
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

    # since no fields are applied the energy offset is simply given by the energies of the kets
    ket_pair_energy_0 = sum(two_atoms_page.ket_config.get_ket_atom(i).get_energy("GHz") for i in range(2))

    two_atoms_page.calculation_config.fast_mode.setChecked(False)
    parameters = ParametersTwoAtoms.from_page(two_atoms_page)
    results = _calculate_two_atoms(parameters)
    energies = np.array(results.energies) + ket_pair_energy_0
    compare_pair_potential_to_reference(energies, np.array(results.ket_overlaps))

    two_atoms_page.calculation_config.fast_mode.setChecked(True)
    parameters = ParametersTwoAtoms.from_page(two_atoms_page)
    results = _calculate_two_atoms(parameters)
    energies = np.array(results.energies) + ket_pair_energy_0
    compare_pair_potential_to_reference(
        energies
    )  # with fast mode, the overlaps are different, so we don't compare them
