import re
from typing import TYPE_CHECKING, Literal, Optional, Union

from PySide6.QtWidgets import (
    QComboBox,
    QLabel,
    QWidget,
)

from pairinteraction_gui.config.base_config import BaseConfig
from pairinteraction_gui.qobjects import (
    DoubleSpinBox,
    HalfIntSpinBox,
    IntSpinBox,
    NamedStackedWidget,
    QnItem,
    WidgetForm,
    WidgetV,
)
from pairinteraction_gui.utils import DatabaseMissingError, NoStateFoundError, get_ket_atom
from pairinteraction_gui.worker import threaded

if TYPE_CHECKING:
    import pairinteraction.real as pi


AVAILABLE_SPECIES = [
    "Rb",
    "Li",
    "Na",
    "K",
    "Cs",
    "Sr88_singlet",
    "Sr88_triplet",
    "Sr87_mqdt",
    "Sr88_mqdt",
    "Yb171_mqdt",
    "Yb173_mqdt",
    "Yb174_mqdt",
]
SpeciesTypes = Literal["sqdt_duplet", "sqdt_singlet", "sqdt_triplet", "mqdt_halfint", "mqdt_int"]


class KetBaseConfig(BaseConfig):
    """Section for configuring the ket of interest."""

    title = "State of Interest"

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
        self.species_combo: list[QComboBox] = []
        self.stacked_qn: list[NamedStackedWidget[QnBase]] = []
        self.ket_label: list[QLabel] = []

    def postSetupWidget(self) -> None:
        self.layout().addStretch(30)

    def setupOneKetAtom(self) -> None:
        """Set up the UI components for a single ket atom."""
        atom = len(self.species_combo)

        # Species selection
        species_widget = WidgetForm(margin=(5, 15, 5, 5), spacing=10)
        species_combo = QComboBox()
        species_combo.addItems(AVAILABLE_SPECIES)  # TODO get available species from pairinteraction
        species_combo.setToolTip("Select the atomic species")
        species_widget.layout().addRow("Species", species_combo)
        species_combo.currentTextChanged.connect(lambda species, atom=atom: self.on_species_changed(species, atom))
        self.layout().addWidget(species_widget)

        # Create stacked widget for different species configurations
        stacked_qn = NamedStackedWidget[QnBase]()
        stacked_qn.addNamedWidget(QnSQDT(s=0.5), "sqdt_duplet")
        stacked_qn.addNamedWidget(QnSQDT(s=0), "sqdt_singlet")
        stacked_qn.addNamedWidget(QnSQDT(s=1), "sqdt_triplet")
        stacked_qn.addNamedWidget(QnMQDT(m_is_int=True), "mqdt_int")
        stacked_qn.addNamedWidget(QnMQDT(m_is_int=False), "mqdt_halfint")

        for _, widget in stacked_qn.items():
            for item in widget.items:
                item.connectAll(lambda atom=atom: self.on_qnitem_changed(atom))
        self.layout().addWidget(stacked_qn)

        # Add a label to display the current ket
        ket_label = QLabel()
        ket_label.setStyleSheet(self._label_style_sheet)
        ket_label.setWordWrap(True)
        self.layout().addWidget(ket_label)

        # Store the widgets for later access
        self.species_combo.append(species_combo)
        self.stacked_qn.append(stacked_qn)
        self.ket_label.append(ket_label)
        self.on_qnitem_changed(atom)

    def get_species(self, atom: int) -> str:
        """Return the selected species of the ... atom."""
        return self.species_combo[atom].currentText()

    def get_species_type(self, atom: int) -> SpeciesTypes:
        """Return the species type based on the species name of the ... atom."""
        species = self.get_species(atom)
        if "mqdt" in species:
            match = re.search(r"\d+", species)
            if match:
                if int(match.group()) % 2 == 0:
                    return "mqdt_int"
                return "mqdt_halfint"
            else:
                raise ValueError(f"Invalid species name: {species}")
        if "singlet" in species:
            return "sqdt_singlet"
        if "triplet" in species:
            return "sqdt_triplet"
        return "sqdt_duplet"

    def get_quantum_numbers(self, atom: int) -> dict[str, float]:
        """Return the quantum numbers of the ... atom."""
        qn_widget = self.stacked_qn[atom].currentWidget()
        qns = {item.label: item.value() for item in qn_widget.items if item.isChecked()}
        return qns

    def get_ket_atom(self, atom: int) -> "pi.KetAtom":
        """Return the ket of interest of the ... atom."""
        species = self.get_species(atom)
        qns = self.get_quantum_numbers(atom)
        return get_ket_atom(species, **qns)

    def on_species_changed(self, species: str, atom: int) -> None:
        """Handle species selection change."""
        species_type = self.get_species_type(atom)
        self.stacked_qn[atom].setCurrentNamedWidget(species_type)
        self.on_qnitem_changed(atom)

    def on_qnitem_changed(self, atom: int) -> None:
        """Update the ket label with current values."""
        try:
            ket = self.get_ket_atom(atom)
            self.ket_label[atom].setText(str(ket))
            self.ket_label[atom].setStyleSheet(self._label_style_sheet)
        except Exception as err:
            if isinstance(err, NoStateFoundError):
                self.ket_label[atom].setText("No ket found. Please select different quantum numbers.")
            elif isinstance(err, DatabaseMissingError):
                self.ket_label[atom].setText("Database required but not downloaded. Please select a different species.")
            else:
                self.ket_label[atom].setText(str(err))
            self.ket_label[atom].setStyleSheet(self._label_style_sheet_error)


class KetAtomConfig(KetBaseConfig):
    def setupWidget(self) -> None:
        super().setupWidget()
        self.setupOneKetAtom()


class KetLifetimesConfig(KetBaseConfig):
    def setupWidget(self) -> None:
        super().setupWidget()
        self.setupOneKetAtom()

        self.layout().addStretch(2)
        spin_box_T = DoubleSpinBox(
            self, vdefault=300, tooltip="Temperature in Kelvin (0K considers only spontaneous decay)"
        )
        self.item_temperature = QnItem(self, "Temperature", spin_box_T, "K")
        self.layout().addWidget(self.item_temperature)
        spin_box_T.valueChanged.connect(lambda value: self.update_lifetime())

        self.layout().addStretch(2)
        # Add a label to display the lifetime
        self.lifetime_label = QLabel()
        self.lifetime_label.setStyleSheet(self._label_style_sheet)
        self.lifetime_label.setWordWrap(True)
        self.layout().addWidget(self.lifetime_label)

        self.update_lifetime()

    def get_temperature(self) -> float:
        if self.item_temperature.isChecked():
            return self.item_temperature.value()
        return 0

    @threaded
    def update_lifetime(self) -> None:
        try:
            ket = self.get_ket_atom(0)
            temperature = self.get_temperature()
            lifetime = ket.get_lifetime(temperature, temperature_unit="K", unit="mus")
            self.lifetime_label.setText(f"Lifetime: {lifetime} μs")
            self.lifetime_label.setStyleSheet(self._label_style_sheet)
        except Exception:
            self.lifetime_label.setText("Ket not found.")
            self.lifetime_label.setStyleSheet(self._label_style_sheet_error)
            self.page.plotwidget.clear()
        else:
            self.page._thread_calculate()

    def on_qnitem_changed(self, atom: int) -> None:
        super().on_qnitem_changed(atom)
        if hasattr(self, "item_temperature"):  # not yet initialized the first time this method is called
            self.update_lifetime()


class KetPairConfig(KetBaseConfig):
    def setupWidget(self) -> None:
        super().setupWidget()

        self.layout().addWidget(QLabel("<b>Atom 1</b>"))
        self.setupOneKetAtom()

        self.layout().addStretch(2)
        self.layout().addWidget(QLabel("<b>Atom 2</b>"))
        self.setupOneKetAtom()

        # set some better defaults
        n_spinbox = self.stacked_qn[1].getNamedWidget("sqdt_duplet").findChild(IntSpinBox, "n")
        n_spinbox.setValue(81)


class QnBase(WidgetV):
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


class QnSQDT(QnBase):
    """Configuration for atoms using SQDT."""

    def __init__(self, parent: Optional[QWidget] = None, s: Union[int, float] = 0.5) -> None:
        self.s = s
        super().__init__(parent)

    def setupWidget(self) -> None:
        spin_boxes = self._spin_boxes = {}

        spin_boxes["n"] = IntSpinBox(self, vmin=1, vdefault=80, tooltip="Principal quantum number n")

        spin_boxes["l"] = IntSpinBox(self, tooltip="Orbital angular momentum l")

        s = self.s
        if s % 1 == 0:
            spin_boxes["j"] = IntSpinBox(self, vmin=int(s), vdefault=int(s), tooltip="Total angular momentum j")
            spin_boxes["m"] = IntSpinBox(self, vmin=-999, vdefault=0, tooltip="Magnetic quantum number m")
        else:
            spin_boxes["j"] = HalfIntSpinBox(self, vmin=s, vdefault=s, tooltip="Total angular momentum j")
            spin_boxes["m"] = HalfIntSpinBox(self, vmin=-999, vdefault=0.5, tooltip="Magnetic quantum number m")


class QnMQDT(QnBase):
    """Configuration for earth alkali atoms using MQDT."""

    default_deactivated: list[str] = ["l_ryd"]

    def __init__(self, parent: Optional[QWidget] = None, m_is_int: bool = False) -> None:
        self.m_is_int = m_is_int
        super().__init__(parent)

    def setupWidget(self) -> None:
        spin_boxes = self._spin_boxes = {}
        spin_boxes["nu"] = DoubleSpinBox(
            self, vmin=1, vdefault=80, vstep=1, tooltip="Effective principal quantum number nu"
        )
        spin_boxes["s"] = DoubleSpinBox(self, vmin=0, vmax=1, vstep=0.1, tooltip="Spin s")
        spin_boxes["j"] = DoubleSpinBox(self, vstep=1, tooltip="Total angular momentum j")

        if self.m_is_int:
            spin_boxes["f"] = IntSpinBox(self, vmin=0, vdefault=0, tooltip="Total angular momentum f")
            spin_boxes["m"] = IntSpinBox(self, vmin=-999, vdefault=0, tooltip="Magnetic quantum number m")
        else:
            spin_boxes["f"] = HalfIntSpinBox(self, vmin=0.5, vdefault=0.5, tooltip="Total angular momentum f")
            spin_boxes["m"] = HalfIntSpinBox(self, vmin=-999, vdefault=0.5, tooltip="Magnetic quantum number m")

        spin_boxes["l_ryd"] = DoubleSpinBox(
            self, vstep=1, tooltip="Orbital angular momentum l_ryd of the Rydberg electron"
        )
