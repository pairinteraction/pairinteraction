# SPDX-FileCopyrightText: 2025 PairInteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later
from __future__ import annotations

from typing import TYPE_CHECKING, Any, Literal, TypedDict

from PySide6.QtCore import Signal
from PySide6.QtWidgets import QComboBox, QLabel

import pairinteraction as pi
from pairinteraction_gui.app import Application
from pairinteraction_gui.config.base_config import BaseConfig
from pairinteraction_gui.qobjects import (
    NamedStackedWidget,
    QnItemDouble,
    QnItemHalfInt,
    QnItemInt,
    WidgetForm,
    WidgetV,
)
from pairinteraction_gui.theme import label_error_theme, label_theme
from pairinteraction_gui.utils import (
    AVAILABLE_SPECIES,
    DatabaseMissingError,
    NoStateFoundError,
    get_custom_error,
    get_species_type,
)
from pairinteraction_gui.worker import MultiThreadWorker

if TYPE_CHECKING:
    from PySide6.QtWidgets import QWidget

    from pairinteraction_gui.page.lifetimes_page import LifetimesPage
    from pairinteraction_gui.qobjects.item import _QnItem

CORE_SPIN_DICT = {"Sr87": 9 / 2, "Sr88": 0, "Yb171": 1 / 2, "Yb173": 5 / 2, "Yb174": 0}


class QuantumNumbers(TypedDict, total=False):
    n: int
    nu: float
    l: float
    s: float
    j: float
    l_ryd: float
    f: float
    m: float


class KetConfig(BaseConfig):
    """Section for configuring the ket of interest."""

    title = "State of Interest"

    signal_species_changed = Signal(int, str)

    def setupWidget(self) -> None:
        self.species_combo_list: list[QComboBox] = []
        self.stacked_qn_list: list[NamedStackedWidget[QnBase]] = []
        self.ket_label_list: list[QLabel] = []

    @property
    def n_atoms(self) -> int:
        """Return the number of atoms configured."""
        return len(self.species_combo_list)

    def setupOneKetAtom(self) -> None:
        """Set up the UI components for a single ket atom."""
        atom = len(self.species_combo_list)

        # Species selection
        species_widget = WidgetForm()
        species_combo = QComboBox()
        species_combo.addItems(AVAILABLE_SPECIES)  # TODO get available species from pairinteraction
        species_combo.setToolTip("Select the atomic species")
        species_widget.layout().addRow("Species", species_combo)
        species_combo.currentTextChanged.connect(
            lambda species, atom=atom: self.signal_species_changed.emit(atom, species)
        )
        self.signal_species_changed.connect(self.on_species_changed)
        self.layout().addWidget(species_widget)

        stacked_qn = NamedStackedWidget[QnBase]()
        self.layout().addWidget(stacked_qn)

        # Add a label to display the current ket
        ket_label = QLabel()
        ket_label.setStyleSheet(label_theme)
        ket_label.setWordWrap(True)
        self.layout().addWidget(ket_label)

        # Store the widgets for later access
        self.species_combo_list.append(species_combo)
        self.stacked_qn_list.append(stacked_qn)
        self.ket_label_list.append(ket_label)
        self.on_qnitem_changed(atom)

    def get_species(self, atom: int = 0) -> str:
        """Return the selected species of the ... atom."""
        return self.species_combo_list[atom].currentText()

    def get_quantum_numbers(self, atom: int = 0) -> QuantumNumbers:
        """Return the quantum numbers of the ... atom."""
        qn_widget = self.stacked_qn_list[atom].currentWidget()
        return {key: item.value() for key, item in qn_widget.items.items() if item.isChecked()}  # type: ignore [return-value]

    def get_ket_atom(self, atom: int, *, ask_download: bool = False) -> pi.KetAtom:
        """Return the ket of interest of the ... atom."""
        species = self.get_species(atom)
        qns = self.get_quantum_numbers(atom)

        try:
            return pi.KetAtom(species, **qns)
        except Exception as err:
            err = get_custom_error(err)
            if ask_download and isinstance(err, DatabaseMissingError):
                Application.signals.ask_download_database.emit(err.species)
            raise err

    def on_species_changed(self, atom: int, species: str) -> None:
        """Handle species selection change."""
        if species not in self.stacked_qn_list[atom]._widgets:
            qn_widget = QnBase.from_species(species, parent=self)
            self.stacked_qn_list[atom].addNamedWidget(qn_widget, species)
            for item in qn_widget.items.values():
                item.connectAll(lambda atom=atom: self.on_qnitem_changed(atom))  # type: ignore [misc]

        self.stacked_qn_list[atom].setCurrentNamedWidget(species)
        self.on_qnitem_changed(atom)

    def on_qnitem_changed(self, atom: int) -> None:
        """Update the ket label with current values."""
        try:
            ket = self.get_ket_atom(atom, ask_download=True)
            self.ket_label_list[atom].setText(str(ket))
            self.ket_label_list[atom].setStyleSheet(label_theme)
        except Exception as err:
            if isinstance(err, NoStateFoundError):
                self.ket_label_list[atom].setText("No ket found. Please select different quantum numbers.")
            elif isinstance(err, DatabaseMissingError):
                self.ket_label_list[atom].setText(
                    "Database required but not downloaded. Please select a different species."
                )
            else:
                self.ket_label_list[atom].setText(str(err))
            self.ket_label_list[atom].setStyleSheet(label_error_theme)


class KetConfigOneAtom(KetConfig):
    def setupWidget(self) -> None:
        super().setupWidget()
        self.setupOneKetAtom()


class KetConfigLifetimes(KetConfig):
    worker_label: MultiThreadWorker | None = None
    worker_plot: MultiThreadWorker | None = None

    page: LifetimesPage

    def setupWidget(self) -> None:
        super().setupWidget()
        self.setupOneKetAtom()
        self.layout().addSpacing(15)

        self.item_temperature = QnItemDouble(
            self,
            "Temperature",
            vdefault=300,
            unit="K",
            tooltip="Temperature in Kelvin (0K considers only spontaneous decay)",
        )
        self.item_temperature.connectAll(self.update_lifetime_label)
        self.layout().addWidget(self.item_temperature)
        self.layout().addSpacing(15)

        # Add a label to display the lifetime
        self.lifetime_label = QLabel()
        self.lifetime_label.setStyleSheet(label_theme)
        self.lifetime_label.setWordWrap(True)
        self.layout().addWidget(self.lifetime_label)

        self.update_lifetime_label()

    def get_temperature(self) -> float:
        return self.item_temperature.value(default=0)

    def update_lifetime_label(self) -> None:
        if self.worker_label and self.worker_label.isRunning():
            self.worker_label.quit()
            self.worker_label.wait()

        def get_lifetime() -> float:
            ket = self.get_ket_atom(0, ask_download=True)
            temperature = self.get_temperature()
            return ket.get_lifetime(temperature, temperature_unit="K", unit="mus")

        def update_result(lifetime: float) -> None:
            self.lifetime_label.setText(f"Lifetime: {lifetime:.3f} μs")
            self.lifetime_label.setStyleSheet(label_theme)

        def update_error(err: Exception) -> None:
            self.lifetime_label.setText("Ket not found.")
            self.lifetime_label.setStyleSheet(label_error_theme)
            self.page.plotwidget.clear()

        self.worker_label = MultiThreadWorker(get_lifetime)
        self.worker_label.signals.result.connect(update_result)
        self.worker_label.signals.error.connect(update_error)
        self.worker_label.start()

        if self.worker_plot and self.worker_plot.isRunning():
            self.worker_plot.quit()
            self.worker_plot.wait()

        self.worker_plot = MultiThreadWorker(self.page.calculate)
        self.worker_plot.signals.result.connect(self.page.update_plot)
        self.worker_plot.start()

    def on_qnitem_changed(self, atom: int) -> None:
        super().on_qnitem_changed(atom)
        if hasattr(self, "item_temperature"):  # not yet initialized the first time this method is called
            self.update_lifetime_label()


class KetConfigTwoAtoms(KetConfig):
    def setupWidget(self) -> None:
        super().setupWidget()

        self.layout().addWidget(QLabel("<b>Atom 1</b>"))
        self.setupOneKetAtom()
        self.layout().addSpacing(15)

        self.layout().addWidget(QLabel("<b>Atom 2</b>"))
        self.setupOneKetAtom()


class QnBase(WidgetV):
    """Base class for quantum number configuration."""

    margin = (10, 0, 10, 0)
    spacing = 5

    items: dict[str, _QnItem[Any]]

    def postSetupWidget(self) -> None:
        for item in self.items.values():
            self.layout().addWidget(item)

    @classmethod
    def from_species(cls, species: str, parent: QWidget | None = None) -> QnBase:
        """Create a quantum number configuration from the species name."""
        species_type = get_species_type(species)
        if species_type == "sqdt_duplet":
            return QnSQDT(parent, s_type="halfint", s=0.5)
        if species_type == "sqdt_singlet":
            return QnSQDT(parent, s_type="int", s=0)
        if species_type == "sqdt_triplet":
            return QnSQDT(parent, s_type="int", s=1)

        element = species.split("_")[0]
        if species_type == "mqdt_halfint":
            i = CORE_SPIN_DICT.get(element, 0.5)
            return QnMQDT(parent, f_type="halfint", i=i)
        if species_type == "mqdt_int":
            i = CORE_SPIN_DICT.get(element, 0)
            return QnMQDT(parent, f_type="int", i=i)

        raise ValueError(f"Unknown species type: {species_type}")


class QnSQDT(QnBase):
    """Configuration for atoms using SQDT."""

    def __init__(
        self,
        parent: QWidget | None = None,
        *,
        s_type: Literal["int", "halfint"],
        s: float,
    ) -> None:
        assert s_type in ("int", "halfint"), "s_type must be int or halfint"
        self.s_type = s_type
        self.s = s
        self.items = {}
        super().__init__(parent)

    def setupWidget(self) -> None:
        self.items["n"] = QnItemInt(self, "n", vmin=1, vdefault=80, tooltip="Principal quantum number n")
        self.items["l"] = QnItemInt(self, "l", vmin=0, tooltip="Orbital angular momentum l")

        s = self.s
        if self.s_type == "int":
            if s != 0:
                self.items["j"] = QnItemInt(self, "j", vmin=int(s), vdefault=int(s), tooltip="Total angular momentum j")
            self.items["m"] = QnItemInt(self, "m", vmin=-999, vmax=999, vdefault=0, tooltip="Magnetic quantum number m")
        else:
            self.items["j"] = QnItemHalfInt(self, "j", vmin=s, vdefault=s, tooltip="Total angular momentum j")
            self.items["m"] = QnItemHalfInt(
                self, "m", vmin=-999.5, vmax=999.5, vdefault=0.5, tooltip="Magnetic quantum number m"
            )


class QnMQDT(QnBase):
    """Configuration for alkaline-earth atoms using MQDT."""

    def __init__(
        self,
        parent: QWidget | None = None,
        *,
        f_type: Literal["int", "halfint"],
        i: float,
    ) -> None:
        assert f_type in ("int", "halfint"), "f_type must be int or halfint"
        self.f_type = f_type
        self.i = i
        self.items = {}
        super().__init__(parent)

    def setupWidget(self) -> None:
        self.items["nu"] = QnItemDouble(
            self, "nu", vmin=1, vdefault=80, vstep=1, tooltip="Effective principal quantum number nu"
        )
        self.items["s"] = QnItemDouble(self, "s", vmin=0, vmax=1, vstep=0.1, tooltip="Spin s")

        if self.i == 0:
            key = "j"
            description = "Total angular momentum j (j=f for I=0)"
        else:
            key = "f"
            description = "Total angular momentum f"
        if self.f_type == "int":
            self.items[key] = QnItemInt(self, key, vmin=0, vdefault=int(self.i), tooltip=description)
            self.items["m"] = QnItemInt(self, "m", vmin=-999, vmax=999, vdefault=0, tooltip="Magnetic quantum number m")
        else:
            self.items[key] = QnItemHalfInt(self, key, vmin=0.5, vdefault=self.i, tooltip=description)
            self.items["m"] = QnItemHalfInt(
                self, "m", vmin=-999.5, vmax=999.5, vdefault=0.5, tooltip="Magnetic quantum number m"
            )

        self.items["l_ryd"] = QnItemDouble(
            self,
            "l_ryd",
            vmin=0,
            vstep=1,
            tooltip="Orbital angular momentum l_ryd of the Rydberg electron",
            checked=False,
        )
