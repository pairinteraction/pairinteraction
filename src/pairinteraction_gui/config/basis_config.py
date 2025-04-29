# SPDX-FileCopyrightText: 2025 Pairinteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

from typing import TYPE_CHECKING, ClassVar, Literal, Union

from PySide6.QtGui import QShowEvent
from PySide6.QtWidgets import (
    QLabel,
)

from pairinteraction import (
    complex as pi_complex,
    real as pi_real,
)
from pairinteraction_gui.config.base_config import BaseConfig
from pairinteraction_gui.qobjects import DoubleSpinBox, HalfIntSpinBox, IntSpinBox, NamedStackedWidget, QnItem, WidgetV
from pairinteraction_gui.utils import DatabaseMissingError, NoStateFoundError
from pairinteraction_gui.worker import Worker

if TYPE_CHECKING:
    from pairinteraction_gui.page import OneAtomPage, TwoAtomsPage


class BasisConfig(BaseConfig):
    """Section for configuring the basis."""

    title = "Basis"
    page: Union["OneAtomPage", "TwoAtomsPage"]

    _label_style_sheet = """
        background-color: #ffffff;
        color: #000000;
        padding: 12px 16px;
        border-radius: 8px;
        margin: 12px 24px 0px 24px;
        border: 1px solid #aaaaaa;
    """

    _label_style_sheet_error = """
        background-color: #fee2e2;
        color: #991b1b;
        padding: 12px 16px;
        border-radius: 8px;
        margin: 12px 24px 0px 24px;
        border: 1px solid #aaaaaa;
    """

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
            for item in widget.items:
                item.connectAll(lambda atom=atom: self.update_basis_label(atom))  # type: ignore [misc]
        self.layout().addWidget(stacked_basis)

        # Add a label to display the current basis
        basis_label = QLabel()
        basis_label.setStyleSheet(self._label_style_sheet)
        basis_label.setWordWrap(True)
        self.layout().addWidget(basis_label)

        # Store the widgets for later access
        self.stacked_basis.append(stacked_basis)
        self.basis_label.append(basis_label)
        self.update_basis_label(atom)

    def update_basis_label(self, atom: int) -> None:
        worker = Worker(self.get_basis, atom)

        def update_result(basis: Union["pi_real.BasisAtom", "pi_complex.BasisAtom"]) -> None:
            self.basis_label[atom].setText(str(basis) + f"\n  ⇒ Basis consists of {basis.number_of_kets} kets")
            self.basis_label[atom].setStyleSheet(self._label_style_sheet)

        worker.signals.result.connect(update_result)

        def update_error(err: Exception) -> None:
            if isinstance(err, NoStateFoundError):
                self.basis_label[atom].setText("Ket of interest wrong quantum numbers, first fix those.")
            elif isinstance(err, DatabaseMissingError):
                self.basis_label[atom].setText(
                    "Database required but not downloaded. Please select a different species."
                )
            else:
                self.basis_label[atom].setText(str(err))
            self.basis_label[atom].setStyleSheet(self._label_style_sheet_error)

        worker.signals.error.connect(update_error)

        worker.start()

    def get_qn_restrictions(self, atom: int) -> dict[str, tuple[float, float]]:
        """Return the quantum number restrictions to construct a BasisAtom."""
        ket = self.page.ket_config.get_ket_atom(atom)
        basis_widget = self.stacked_basis[atom].currentWidget()
        delta_qns: dict[str, float] = {item.label: item.value() for item in basis_widget.items if item.isChecked()}

        qns: dict[str, tuple[float, float]] = {}
        for key, value in delta_qns.items():
            if value < 0:
                continue
            key = key.replace("Δ", "")
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
        delta_qns: dict[str, float] = {}
        for item in self.stacked_basis[atom].currentWidget().items:
            if item.isChecked():
                delta_qns[item.label.replace("Δ", "")] = item.value()
        return delta_qns

    def showEvent(self, event: QShowEvent) -> None:
        """Update the basis label when the widget is shown."""
        super().showEvent(event)
        ket_of_interest = self.page.ket_config
        for atom in range(len(ket_of_interest.species_combo)):
            species_type = ket_of_interest.get_species_type(atom)
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
        self.delta_pair_energy = DoubleSpinBox(
            self, vmin=0, vdefault=5, tooltip="Restriction for the pair energy difference to the state of interest"
        )
        self.layout().addWidget(QnItem(self, "ΔEnergy", self.delta_pair_energy, "GHz"))

        self.basis_pair_label = QLabel()
        self.basis_pair_label.setStyleSheet(self._label_style_sheet)
        self.basis_pair_label.setWordWrap(True)
        self.layout().addWidget(self.basis_pair_label)

    def update_basis_pair_label(self, basis_pair_label: str) -> None:
        """Update the quantum state label with current values."""
        self.basis_pair_label.setText(basis_pair_label)
        self.basis_pair_label.setStyleSheet(self._label_style_sheet)

    def clear_basis_pair_label(self) -> None:
        """Clear the basis pair label."""
        self.basis_pair_label.setText("")


class RestrictionsBase(WidgetV):
    """Base class for quantum number configuration."""

    margin = (10, 0, 10, 0)
    spacing = 5

    default_deactivated: ClassVar[list[str]] = []
    _spin_boxes: dict[str, Union[IntSpinBox, HalfIntSpinBox, DoubleSpinBox]]
    items: list[QnItem]

    def postSetupWidget(self) -> None:
        self.items = []
        for key, spin_box in self._spin_boxes.items():
            unit = "GHz" if "Energy" in key else ""
            checked = key not in self.default_deactivated
            item = QnItem(self, key, spin_box, unit, checked)
            self.items.append(item)
            self.layout().addWidget(item)


class RestrictionsSQDT(RestrictionsBase):
    """Configuration for alkali atoms using SQDT."""

    default_deactivated: ClassVar[list[str]] = ["Δj", "Δm"]

    def setupWidget(self) -> None:
        spin_boxes = self._spin_boxes = {}
        spin_boxes["Δn"] = IntSpinBox(self, vdefault=3, tooltip="Restriction for the Principal quantum number n")
        spin_boxes["Δl"] = IntSpinBox(self, vdefault=2, tooltip="Restriction for the Orbital angular momentum l")
        spin_boxes["Δj"] = IntSpinBox(self, tooltip="Restriction for the Total angular momentum j")
        spin_boxes["Δm"] = IntSpinBox(self, tooltip="Restriction for the Magnetic quantum number m")


class RestrictionsMQDT(RestrictionsBase):
    """Configuration for alkali atoms using SQDT."""

    default_deactivated: ClassVar[list[str]] = ["Δf", "Δm", "Δl_ryd"]

    def setupWidget(self) -> None:
        spin_boxes = self._spin_boxes = {}
        spin_boxes["Δnu"] = DoubleSpinBox(
            self, vdefault=4, tooltip="Restriction for the Effective principal quantum number nu"
        )
        spin_boxes["Δs"] = DoubleSpinBox(self, vdefault=0.5, tooltip="Restriction for the Spin s")
        spin_boxes["Δj"] = DoubleSpinBox(self, vdefault=3, tooltip="Restriction for the Total angular momentum j")

        spin_boxes["Δf"] = IntSpinBox(self, vdefault=5, tooltip="Restriction for the Total angular momentum f")
        spin_boxes["Δm"] = IntSpinBox(self, vdefault=5, tooltip="Restriction for the Magnetic quantum number m")

        # TODO this is not yet implemented
        # spin_boxes["Δl_ryd"] = DoubleSpinBox(
        #     self, vdefault=3, tooltip="Restriction for the Orbital angular momentum l_ryd of the Rydberg electron"
        # )
