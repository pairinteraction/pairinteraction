# SPDX-FileCopyrightText: 2025 PairInteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

from typing import TYPE_CHECKING, Any, Literal, Optional, TypedDict, Union

from PySide6.QtGui import QShowEvent
from PySide6.QtWidgets import (
    QLabel,
    QWidget,
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
    from pairinteraction_gui.config.ket_config import QuantumNumbers
    from pairinteraction_gui.page import OneAtomPage, TwoAtomsPage
    from pairinteraction_gui.qobjects.item import _QnItem


class QuantumNumberRestrictions(TypedDict, total=False):
    n: tuple[int, int]
    nu: tuple[float, float]
    l: tuple[float, float]
    s: tuple[float, float]
    j: tuple[float, float]
    l_ryd: tuple[float, float]
    f: tuple[float, float]
    m: tuple[float, float]


class BasisConfig(BaseConfig):
    """Section for configuring the basis."""

    title = "Basis"
    page: Union["OneAtomPage", "TwoAtomsPage"]

    def setupWidget(self) -> None:
        self.stacked_basis_list: list[NamedStackedWidget[RestrictionsBase]] = []
        self.basis_label_list: list[QLabel] = []

    def setupOneBasisAtom(self) -> None:
        """Set up the UI components for a single basis atom."""
        atom = len(self.stacked_basis_list)

        stacked_basis = NamedStackedWidget[RestrictionsBase]()
        self.layout().addWidget(stacked_basis)

        # Add a label to display the current basis
        basis_label = QLabel()
        basis_label.setStyleSheet(label_theme)
        basis_label.setWordWrap(True)
        self.layout().addWidget(basis_label)

        # Store the widgets for later access
        self.stacked_basis_list.append(stacked_basis)
        self.basis_label_list.append(basis_label)
        self.update_basis_label(atom)

    def update_basis_label(self, atom: int) -> None:
        worker = MultiThreadWorker(self.get_basis, atom)

        def update_result(basis: Union["pi_real.BasisAtom", "pi_complex.BasisAtom"]) -> None:
            self.basis_label_list[atom].setText(str(basis) + f"\n  ⇒ Basis consists of {basis.number_of_kets} kets")
            self.basis_label_list[atom].setStyleSheet(label_theme)

        worker.signals.result.connect(update_result)

        def update_error(err: Exception) -> None:
            if isinstance(err, NoStateFoundError):
                self.basis_label_list[atom].setText("Ket of interest wrong quantum numbers, first fix those.")
            elif isinstance(err, DatabaseMissingError):
                self.basis_label_list[atom].setText(
                    "Database required but not downloaded. Please select a different state of interest."
                )
            else:
                self.basis_label_list[atom].setText(str(err))
            self.basis_label_list[atom].setStyleSheet(label_error_theme)

        worker.signals.error.connect(update_error)

        worker.start()

    def get_quantum_number_restrictions(self, atom: int) -> QuantumNumberRestrictions:
        """Return the quantum number restrictions to construct a BasisAtom."""
        ket = self.page.ket_config.get_ket_atom(atom)
        qns = self.page.ket_config.get_quantum_numbers(atom)
        delta_qns = self.get_quantum_number_deltas(atom)

        qn_restrictions: dict[str, tuple[float, float]] = {}
        for key, delta in delta_qns.items():
            if key in qns:
                qn = qns[key]  # type: ignore [literal-required]
            elif hasattr(ket, key):
                qn = getattr(ket, key)
            else:
                raise ValueError(f"Quantum number {key} not found in quantum_numbers or KetAtom.")
            qn_restrictions[key] = (qn - delta, qn + delta)

        return qn_restrictions  # type: ignore [return-value]

    def get_basis(
        self, atom: int, dtype: Literal["real", "complex"] = "real"
    ) -> Union["pi_real.BasisAtom", "pi_complex.BasisAtom"]:
        """Return the basis of interest."""
        ket = self.page.ket_config.get_ket_atom(atom)
        qn_restrictions = self.get_quantum_number_restrictions(atom)
        if dtype == "real":
            return pi_real.BasisAtom(ket.species, **qn_restrictions)
        return pi_complex.BasisAtom(ket.species, **qn_restrictions)

    def get_quantum_number_deltas(self, atom: int = 0) -> "QuantumNumbers":
        """Return the quantum number deltas for the basis of interest."""
        basis_widget = self.stacked_basis_list[atom].currentWidget()
        return {key: item.value() for key, item in basis_widget.items.items() if item.isChecked()}  # type: ignore [return-value]

    def on_species_changed(self, atom: int, species: str) -> None:
        """Handle species selection change."""
        if species not in self.stacked_basis_list[atom]._widgets:
            restrictions_widget = RestrictionsBase.from_species(species, parent=self)
            self.stacked_basis_list[atom].addNamedWidget(restrictions_widget, species)
            for _, item in restrictions_widget.items.items():
                item.connectAll(lambda atom=atom: self.update_basis_label(atom))  # type: ignore [misc]

        self.stacked_basis_list[atom].setCurrentNamedWidget(species)
        self.update_basis_label(atom)

    def showEvent(self, event: QShowEvent) -> None:
        super().showEvent(event)
        for i in range(len(self.stacked_basis_list)):
            self.update_basis_label(i)


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

    @classmethod
    def from_species(cls, species: str, parent: Optional[QWidget] = None) -> "RestrictionsBase":
        """Create a quantum number restriction configuration from the species name."""
        species_type = get_species_type(species)
        if species_type == "sqdt_duplet":
            return RestrictionsSQDT(parent, s_type="halfint", s=0.5)
        if species_type == "sqdt_singlet":
            return RestrictionsSQDT(parent, s_type="int", s=0)
        if species_type == "sqdt_triplet":
            return RestrictionsSQDT(parent, s_type="int", s=1)
        if species_type == "mqdt_halfint":
            return RestrictionsMQDT(parent, f_type="halfint", i=0.5)
        if species_type == "mqdt_int":
            return RestrictionsMQDT(parent, f_type="int", i=0)

        raise ValueError(f"Unknown species type: {species_type}")


class RestrictionsSQDT(RestrictionsBase):
    """Configuration atoms using SQDT."""

    def __init__(
        self,
        parent: Optional[QWidget] = None,
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
        self.items["n"] = QnItemInt(self, "Δn", vdefault=3, tooltip="Restriction for the principal quantum number n")
        self.items["l"] = QnItemInt(self, "Δl", vdefault=2, tooltip="Restriction for the orbital angular momentum l")
        if self.s != 0:
            self.items["j"] = QnItemInt(
                self, "Δj", tooltip="Restriction for the total angular momentum j", checked=False
            )
        self.items["m"] = QnItemInt(self, "Δm", tooltip="Restriction for the magnetic quantum number m", checked=False)


class RestrictionsMQDT(RestrictionsBase):
    """Configuration for alkali atoms using SQDT."""

    def __init__(
        self,
        parent: Optional[QWidget] = None,
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
            self, "Δnu", vdefault=4, tooltip="Restriction for the effective principal quantum number nu"
        )
        self.items["s"] = QnItemDouble(self, "Δs", vdefault=0.5, tooltip="Restriction for the spin s")

        if self.i == 0:
            key = "j"
            description = "Restriction for the  total angular momentum j (j=f for I=0)"
        else:
            key = "f"
            description = "Restriction for the  total angular momentum f"
        self.items[key] = QnItemDouble(self, "Δ" + key, vdefault=3, tooltip=description)

        self.items["m"] = QnItemInt(
            self, "Δm", vdefault=5, tooltip="Restriction for the magnetic quantum number m", checked=False
        )

        self.items["l_ryd"] = QnItemDouble(
            self,
            "Δl_ryd",
            vdefault=3,
            tooltip="Restriction for the orbital angular momentum l_ryd of the Rydberg electron",
            checked=False,
        )
