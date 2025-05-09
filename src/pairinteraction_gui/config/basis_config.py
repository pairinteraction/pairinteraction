# SPDX-FileCopyrightText: 2025 Pairinteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

from typing import TYPE_CHECKING, Any, Literal, Union

from PySide6.QtWidgets import (
    QLabel,
)

from pairinteraction import (
    complex as pi_complex,
    real as pi_real,
)
from pairinteraction_gui.config.base_config import BaseConfig
from pairinteraction_gui.qobjects import NamedStackedWidget, QnItemDouble, QnItemInt, WidgetV
from pairinteraction_gui.qobjects.item import RangeItem
from pairinteraction_gui.theme import label_error_theme, label_theme
from pairinteraction_gui.utils import DatabaseMissingError, NoStateFoundError, get_species_type
from pairinteraction_gui.worker import MultiThreadWorker

if TYPE_CHECKING:
    from pairinteraction_gui.page import OneAtomPage, TwoAtomsPage
    from pairinteraction_gui.qobjects.item import _QnItem


class BasisConfig(BaseConfig):
    """Section for configuring the basis."""

    title = "Basis"
    page: Union["OneAtomPage", "TwoAtomsPage"]

    def setupWidget(self) -> None:
        self.stacked_basis: list[NamedStackedWidget[RestrictionsBase]] = []
        self.basis_label: list[QLabel] = []

    def setupOneBasisAtom(self) -> None:
        """Set up the UI components for a single basis atom."""
        atom = len(self.stacked_basis)

        # Create stacked widget for different species configurations
        stacked_basis = NamedStackedWidget[RestrictionsBase]()
        stacked_basis.addNamedWidget(RestrictionsSQDT(), "sqdt")
        stacked_basis.addNamedWidget(RestrictionsMQDT(), "mqdt")

        for _, widget in stacked_basis.items():
            for _, item in widget.items.items():
                item.connectAll(lambda atom=atom: self.update_basis_label(atom))  # type: ignore [misc]
        self.layout().addWidget(stacked_basis)

        # Add a label to display the current basis
        basis_label = QLabel()
        basis_label.setStyleSheet(label_theme)
        basis_label.setWordWrap(True)
        self.layout().addWidget(basis_label)

        # Store the widgets for later access
        self.stacked_basis.append(stacked_basis)
        self.basis_label.append(basis_label)
        self.update_basis_label(atom)

    def update_basis_label(self, atom: int) -> None:
        worker = MultiThreadWorker(self.get_basis, atom)

        def update_result(basis: Union["pi_real.BasisAtom", "pi_complex.BasisAtom"]) -> None:
            self.basis_label[atom].setText(str(basis) + f"\n  ⇒ Basis consists of {basis.number_of_kets} kets")
            self.basis_label[atom].setStyleSheet(label_theme)

        worker.signals.result.connect(update_result)

        def update_error(err: Exception) -> None:
            if isinstance(err, NoStateFoundError):
                self.basis_label[atom].setText("Ket of interest wrong quantum numbers, first fix those.")
            elif isinstance(err, DatabaseMissingError):
                self.basis_label[atom].setText(
                    "Database required but not downloaded. Please select a different state of interest."
                )
            else:
                self.basis_label[atom].setText(str(err))
            self.basis_label[atom].setStyleSheet(label_error_theme)

        worker.signals.error.connect(update_error)

        worker.start()

    def get_qn_restrictions(self, atom: int) -> dict[str, tuple[float, float]]:
        """Return the quantum number restrictions to construct a BasisAtom."""
        ket = self.page.ket_config.get_ket_atom(atom)
        basis_widget = self.stacked_basis[atom].currentWidget()
        delta_qns: dict[str, float] = {
            key: item.value() for key, item in basis_widget.items.items() if item.isChecked()
        }

        qns: dict[str, tuple[float, float]] = {}
        for key, value in delta_qns.items():
            qn: float = getattr(ket, key)
            qns[key] = (qn - value, qn + value)

        return qns

    def get_basis(
        self, atom: int, dtype: Literal["real", "complex"] = "real"
    ) -> Union["pi_real.BasisAtom", "pi_complex.BasisAtom"]:
        """Return the basis of interest."""
        ket = self.page.ket_config.get_ket_atom(atom)
        qn_restrictions = self.get_qn_restrictions(atom)
        if dtype == "real":
            return pi_real.BasisAtom(ket.species, **qn_restrictions)  # type: ignore [arg-type]
        return pi_complex.BasisAtom(ket.species, **qn_restrictions)  # type: ignore [arg-type]

    def get_quantum_number_deltas(self, atom: int = 0) -> dict[str, float]:
        """Return the quantum number deltas for the basis of interest."""
        stacked_basis = self.stacked_basis[atom].currentWidget()
        return {key: item.value() for key, item in stacked_basis.items.items() if item.isChecked()}

    def on_species_changed(self, atom: int, species: str) -> None:
        """Handle species selection change."""
        species_type = get_species_type(species)
        if "mqdt" in species_type:
            self.stacked_basis[atom].setCurrentNamedWidget("mqdt")
        else:
            self.stacked_basis[atom].setCurrentNamedWidget("sqdt")
        self.update_basis_label(atom)


class BasisConfigOneAtom(BasisConfig):
    def setupWidget(self) -> None:
        super().setupWidget()
        self.setupOneBasisAtom()


class BasisConfigTwoAtoms(BasisConfig):
    def setupWidget(self) -> None:
        super().setupWidget()

        self.layout().addWidget(QLabel("<b>Atom 1</b>"))
        self.setupOneBasisAtom()
        self.layout().addSpacing(15)

        self.layout().addWidget(QLabel("<b>Atom 2</b>"))
        self.setupOneBasisAtom()
        self.layout().addSpacing(15)

        self.layout().addWidget(QLabel("<b>Pair Basis Restrictions</b>"))
        self.pair_delta_energy = QnItemDouble(
            self,
            "ΔEnergy",
            vdefault=5,
            vmin=0,
            unit="GHz",
            tooltip="Restriction for the pair energy difference to the state of interest",
        )
        self.layout().addWidget(self.pair_delta_energy)
        self.pair_m_range = RangeItem(
            self,
            "Total m",
            tooltip_label="pair total angular momentum m",
            checked=False,
        )
        self.layout().addWidget(self.pair_m_range)

        self.basis_pair_label = QLabel()
        self.basis_pair_label.setStyleSheet(label_theme)
        self.basis_pair_label.setWordWrap(True)
        self.layout().addWidget(self.basis_pair_label)

    def update_basis_pair_label(self, basis_pair_label: str) -> None:
        """Update the quantum state label with current values."""
        self.basis_pair_label.setText(basis_pair_label)
        self.basis_pair_label.setStyleSheet(label_theme)

    def clear_basis_pair_label(self) -> None:
        """Clear the basis pair label."""
        self.basis_pair_label.setText("")


class RestrictionsBase(WidgetV):
    """Base class for quantum number configuration."""

    margin = (10, 0, 10, 0)
    spacing = 5

    items: dict[str, "_QnItem[Any]"]

    def postSetupWidget(self) -> None:
        for _key, item in self.items.items():
            self.layout().addWidget(item)


class RestrictionsSQDT(RestrictionsBase):
    """Configuration for alkali atoms using SQDT."""

    def setupWidget(self) -> None:
        self.items = {}
        self.items["n"] = QnItemInt(self, "Δn", vdefault=3, tooltip="Restriction for the Principal quantum number n")
        self.items["l"] = QnItemInt(self, "Δl", vdefault=2, tooltip="Restriction for the Orbital angular momentum l")
        self.items["j"] = QnItemInt(self, "Δj", tooltip="Restriction for the Total angular momentum j", checked=False)
        self.items["m"] = QnItemInt(self, "Δm", tooltip="Restriction for the Magnetic quantum number m", checked=False)


class RestrictionsMQDT(RestrictionsBase):
    """Configuration for alkali atoms using SQDT."""

    def setupWidget(self) -> None:
        self.items = {}
        self.items["nu"] = QnItemDouble(
            self, "Δnu", vdefault=4, tooltip="Restriction for the Effective principal quantum number nu"
        )
        self.items["s"] = QnItemDouble(self, "Δs", vdefault=0.5, tooltip="Restriction for the Spin s")
        self.items["j"] = QnItemDouble(self, "Δj", vdefault=3, tooltip="Restriction for the Total angular momentum j")

        self.items["f"] = QnItemInt(
            self, "Δf", vdefault=5, tooltip="Restriction for the Total angular momentum f", checked=False
        )
        self.items["m"] = QnItemInt(
            self, "Δm", vdefault=5, tooltip="Restriction for the Magnetic quantum number m", checked=False
        )
