# SPDX-FileCopyrightText: 2025 Pairinteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

from typing import TYPE_CHECKING, Any, Optional

from PySide6.QtCore import Signal
from PySide6.QtWidgets import QComboBox, QLabel, QWidget

import pairinteraction.real as pi
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
    from pairinteraction_gui.page.lifetimes_page import LifetimesPage
    from pairinteraction_gui.qobjects.item import _QnItem


class KetConfig(BaseConfig):
    """Section for configuring the ket of interest."""

    title = "State of Interest"

    signal_species_changed = Signal(int, str)

    def setupWidget(self) -> None:
        self.species_combo: list[QComboBox] = []
        self.stacked_qn: list[NamedStackedWidget[QnBase]] = []
        self.ket_label: list[QLabel] = []

    @property
    def n_atoms(self) -> int:
        """Return the number of atoms configured."""
        return len(self.species_combo)

    def setupOneKetAtom(self) -> None:
        """Set up the UI components for a single ket atom."""
        atom = len(self.species_combo)

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

        # Create stacked widget for different species configurations
        stacked_qn = NamedStackedWidget[QnBase]()
        stacked_qn.addNamedWidget(QnSQDT(s=0.5), "sqdt_duplet")
        stacked_qn.addNamedWidget(QnSQDT(s=0), "sqdt_singlet")
        stacked_qn.addNamedWidget(QnSQDT(s=1), "sqdt_triplet")
        stacked_qn.addNamedWidget(QnMQDT(m_is_int=True), "mqdt_int")
        stacked_qn.addNamedWidget(QnMQDT(m_is_int=False), "mqdt_halfint")

        for _, widget in stacked_qn.items():
            for item in widget.items.values():
                item.connectAll(lambda atom=atom: self.on_qnitem_changed(atom))  # type: ignore [misc]
        self.layout().addWidget(stacked_qn)

        # Add a label to display the current ket
        ket_label = QLabel()
        ket_label.setStyleSheet(label_theme)
        ket_label.setWordWrap(True)
        self.layout().addWidget(ket_label)

        # Store the widgets for later access
        self.species_combo.append(species_combo)
        self.stacked_qn.append(stacked_qn)
        self.ket_label.append(ket_label)
        self.on_qnitem_changed(atom)

    def get_species(self, atom: int = 0) -> str:
        """Return the selected species of the ... atom."""
        return self.species_combo[atom].currentText()

    def get_quantum_numbers(self, atom: int = 0) -> dict[str, float]:
        """Return the quantum numbers of the ... atom."""
        qn_widget = self.stacked_qn[atom].currentWidget()
        return {key: item.value() for key, item in qn_widget.items.items() if item.isChecked()}

    def get_ket_atom(self, atom: int, *, ask_download: bool = False) -> "pi.KetAtom":
        """Return the ket of interest of the ... atom."""
        species = self.get_species(atom)
        qns = self.get_quantum_numbers(atom)

        try:
            return pi.KetAtom(species, **qns)  # type: ignore [arg-type]
        except Exception as err:
            err = get_custom_error(err)
            if ask_download and isinstance(err, DatabaseMissingError):
                Application.signals.ask_download_database.emit(err.species)
            raise err

    def on_species_changed(self, atom: int, species: str) -> None:
        """Handle species selection change."""
        species_type = get_species_type(species)
        self.stacked_qn[atom].setCurrentNamedWidget(species_type)
        self.on_qnitem_changed(atom)

    def on_qnitem_changed(self, atom: int) -> None:
        """Update the ket label with current values."""
        try:
            ket = self.get_ket_atom(atom, ask_download=True)
            self.ket_label[atom].setText(str(ket))
            self.ket_label[atom].setStyleSheet(label_theme)
        except Exception as err:
            if isinstance(err, NoStateFoundError):
                self.ket_label[atom].setText("No ket found. Please select different quantum numbers.")
            elif isinstance(err, DatabaseMissingError):
                self.ket_label[atom].setText("Database required but not downloaded. Please select a different species.")
            else:
                self.ket_label[atom].setText(str(err))
            self.ket_label[atom].setStyleSheet(label_error_theme)


class KetConfigOneAtom(KetConfig):
    def setupWidget(self) -> None:
        super().setupWidget()
        self.setupOneKetAtom()


class KetConfigLifetimes(KetConfig):
    worker_label: Optional[MultiThreadWorker] = None
    worker_plot: Optional[MultiThreadWorker] = None

    page: "LifetimesPage"

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

    items: dict[str, "_QnItem[Any]"]

    def postSetupWidget(self) -> None:
        for item in self.items.values():
            self.layout().addWidget(item)


class QnSQDT(QnBase):
    """Configuration for atoms using SQDT."""

    def __init__(self, parent: Optional[QWidget] = None, s: float = 0.5) -> None:
        self.s = s
        self.items = {}
        super().__init__(parent)

    def setupWidget(self) -> None:
        self.items["n"] = QnItemInt(self, "n", vmin=1, vdefault=80, tooltip="Principal quantum number n")
        self.items["l"] = QnItemInt(self, "l", vmin=0, tooltip="Orbital angular momentum l")

        s = self.s
        if s % 1 == 0:
            self.items["j"] = QnItemInt(self, "j", vmin=int(s), vdefault=int(s), tooltip="Total angular momentum j")
            self.items["m"] = QnItemInt(self, "m", vmin=-999, vmax=999, vdefault=0, tooltip="Magnetic quantum number m")
        else:
            self.items["j"] = QnItemHalfInt(self, "j", vmin=s, vdefault=s, tooltip="Total angular momentum j")
            self.items["m"] = QnItemHalfInt(
                self, "m", vmin=-999.5, vmax=999.5, vdefault=0.5, tooltip="Magnetic quantum number m"
            )


class QnMQDT(QnBase):
    """Configuration for earth alkali atoms using MQDT."""

    def __init__(self, parent: Optional[QWidget] = None, m_is_int: bool = False) -> None:
        self.m_is_int = m_is_int
        self.items = {}
        super().__init__(parent)

    def setupWidget(self) -> None:
        self.items["nu"] = QnItemDouble(
            self, "nu", vmin=1, vdefault=80, vstep=1, tooltip="Effective principal quantum number nu"
        )
        self.items["s"] = QnItemDouble(self, "s", vmin=0, vmax=1, vstep=0.1, tooltip="Spin s")
        self.items["j"] = QnItemDouble(self, "j", vmin=0, vstep=1, tooltip="Total angular momentum j")

        if self.m_is_int:
            self.items["f"] = QnItemInt(self, "f", vmin=0, vdefault=0, tooltip="Total angular momentum f")
            self.items["m"] = QnItemInt(self, "m", vmin=-999, vmax=999, vdefault=0, tooltip="Magnetic quantum number m")
        else:
            self.items["f"] = QnItemHalfInt(self, "f", vmin=0.5, vdefault=0.5, tooltip="Total angular momentum f")
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
