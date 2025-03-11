from typing import TYPE_CHECKING, Literal, Union

from PySide6.QtGui import QShowEvent
from PySide6.QtWidgets import (
    QLabel,
)

from pairinteraction_gui.config.base_config import BaseConfig
from pairinteraction_gui.qobjects import DoubleSpinBox, HalfIntSpinBox, IntSpinBox, NamedStackedWidget, QnItem, WidgetV
from pairinteraction_gui.utils import DatabaseMissingError, NoStateFoundError, get_basis_atom
from pairinteraction_gui.worker import threaded

if TYPE_CHECKING:
    from pairinteraction import (
        complex as pi_complex,
        real as pi_real,
    )
    from pairinteraction._wrapped.basis.BasisPair import BasisPairBase
    from pairinteraction_gui.page import SystemAtomPage, SystemPairPage


class BasisBaseConfig(BaseConfig):
    """Section for configuring the basis."""

    title = "Basis"
    page: Union["SystemAtomPage", "SystemPairPage"]

    _label_style_sheet = """
        background-color: #f8f9fa;
        color: #212529;
        padding: 12px 16px;
        border-radius: 8px;
        margin: 12px 24px 0px 24px;
        border: 1px solid #212529;
    """

    _label_style_sheet_error = """
        background-color: #fee2e2;
        color: #991b1b;
        padding: 12px 16px;
        border-radius: 8px;
        margin: 12px 24px 0px 24px;
        border: 1px solid #dc2626;
    """

    def setupWidget(self) -> None:
        self.stacked_basis: list[NamedStackedWidget[RestrictionsBase]] = []
        self.basis_label: list[QLabel] = []

    def postSetupWidget(self) -> None:
        self.layout().addStretch(30)

    def setupOneBasisAtom(self) -> None:
        """Set up the UI components for a single basis atom."""
        atom = len(self.stacked_basis)

        # Create stacked widget for different species configurations
        stacked_basis = NamedStackedWidget[RestrictionsBase]()
        stacked_basis.addNamedWidget(RestrictionsSQDT(), "sqdt")
        stacked_basis.addNamedWidget(RestrictionsMQDT(), "mqdt")

        for _, widget in stacked_basis.items():
            for item in widget.items:
                item.connectAll(lambda atom=atom: self.update_basis_label(atom))
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

    @threaded
    def update_basis_label(self, atom: int) -> None:
        """Update the quantum state label with current values."""
        try:
            basis = self.get_basis(atom)
            self.basis_label[atom].setText(str(basis) + f"\n  ⇒ Basis consists of {basis.number_of_kets} kets")
            self.basis_label[atom].setStyleSheet(self._label_style_sheet)
        except Exception as err:
            if isinstance(err, NoStateFoundError):
                self.basis_label[atom].setText("Ket of interest wrong quantum numbers, first fix those.")
            elif isinstance(err, DatabaseMissingError):
                self.basis_label[atom].setText(
                    "Database required but not downloaded. Please select a different species."
                )
            else:
                self.basis_label[atom].setText(str(err))
            self.basis_label[atom].setStyleSheet(self._label_style_sheet_error)

    def get_basis(
        self, atom: int, dtype: Literal["real", "complex"] = "real"
    ) -> Union["pi_real.BasisAtom", "pi_complex.BasisAtom"]:
        """Return the basis of interest."""
        basis_widget = self.stacked_basis[atom].currentWidget()
        delta_qns: dict[str, Union[float, int]] = {
            item.label: item.value() for item in basis_widget.items if item.isChecked()
        }
        ket = self.page.ket_config.get_ket_atom(atom)
        return get_basis_atom(ket, **delta_qns, dtype=dtype)

    def showEvent(self, event: "QShowEvent") -> None:
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


class BasisAtomConfig(BasisBaseConfig):
    def setupWidget(self) -> None:
        super().setupWidget()
        self.setupOneBasisAtom()


class BasisPairConfig(BasisBaseConfig):
    def setupWidget(self) -> None:
        super().setupWidget()

        self.layout().addWidget(QLabel("<b>Atom 1</b>"))
        self.setupOneBasisAtom()

        self.layout().addStretch(5)
        self.layout().addWidget(QLabel("<b>Atom 2</b>"))
        self.setupOneBasisAtom()

        self.layout().addStretch(5)
        self.layout().addWidget(QLabel("<b>Pair Basis Restrictions</b>"))
        self.delta_pair_energy = DoubleSpinBox(
            self, vmin=0, vdefault=3, tooltip="Restriction for the pair energy difference to the state of interest"
        )
        self.layout().addWidget(QnItem(self, "ΔEnergy", self.delta_pair_energy, "GHz"))

        self.basis_pair_label = QLabel()
        self.basis_pair_label.setStyleSheet(self._label_style_sheet)
        self.basis_pair_label.setWordWrap(True)
        self.layout().addWidget(self.basis_pair_label)

    def update_basis_pair_label(self, basis_pair: "BasisPairBase") -> None:
        """Update the quantum state label with current values."""
        self.basis_pair_label.setText(str(basis_pair) + f"\n  ⇒ Basis consists of {basis_pair.number_of_kets} kets")
        self.basis_pair_label.setStyleSheet(self._label_style_sheet)

    def clear_basis_pair_label(self) -> None:
        """Clear the basis pair label."""
        self.basis_pair_label.setText("")


class RestrictionsBase(WidgetV):
    """Base class for quantum number configuration."""

    margin = (10, 0, 10, 0)
    spacing = 5

    default_deactivated: list[str] = []
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

    default_deactivated = ["Δj", "Δm"]

    def setupWidget(self) -> None:
        spin_boxes = self._spin_boxes = {}
        spin_boxes["Δn"] = IntSpinBox(self, vdefault=3, tooltip="Restriction for the Principal quantum number n")
        spin_boxes["Δl"] = IntSpinBox(self, vdefault=2, tooltip="Restriction for the Orbital angular momentum l")
        spin_boxes["Δj"] = IntSpinBox(self, tooltip="Restriction for the Total angular momentum j")
        spin_boxes["Δm"] = IntSpinBox(self, tooltip="Restriction for the Magnetic quantum number m")


class RestrictionsMQDT(RestrictionsBase):
    """Configuration for alkali atoms using SQDT."""

    default_deactivated = ["Δf", "Δm", "Δl_ryd"]

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
