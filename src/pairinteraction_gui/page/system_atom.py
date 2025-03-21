import logging
from collections.abc import Sequence
from pathlib import Path
from typing import TYPE_CHECKING, Any

import nbformat
from nbconvert import PythonExporter
from PySide6.QtWidgets import QFileDialog

import pairinteraction._wrapped as pi
from pairinteraction_gui.config import BasisAtomConfig, KetAtomConfig, SystemAtomConfig
from pairinteraction_gui.page.base_page import SimulationPage
from pairinteraction_gui.plotwidget.plotwidget import PlotEnergies

if TYPE_CHECKING:
    from numpy.typing import NDArray

logger = logging.getLogger(__name__)


class SystemAtomPage(SimulationPage):
    """Page for configuring and analyzing single atom systems."""

    title = "One Atom"
    tooltip = "Configure and analyze single atom systems"

    plotwidget: PlotEnergies

    def setupWidget(self) -> None:
        self.plotwidget = PlotEnergies(self)
        self.layout().addWidget(self.plotwidget)
        super().setupWidget()

        # all attributes of instance BaseConfig will be added to the toolbox in postSetupWidget
        self.ket_config = KetAtomConfig(self)
        self.basis_config = BasisAtomConfig(self)
        self.system_atom_config = SystemAtomConfig(self)

    def calculate(self) -> None:
        super().calculate()

        self.fields = self.system_atom_config.get_fields()
        self.systems = self.system_atom_config.get_systems(0)
        self.ket = self.ket_config.get_ket_atom(0)
        kwargs: dict[str, Any] = {}
        if self.plotwidget.energy_range.isChecked():
            _energies = self.plotwidget.energy_range.values()
            kwargs["energy_range"] = (self.ket.get_energy("GHz") + v for v in _energies)
            kwargs["energy_unit"] = "GHz"
        if self.plotwidget.fast_mode.isChecked():
            kwargs["diagonalizer"] = "lapacke_evr"
            kwargs["float_type"] = "float32"

        pi.diagonalize(self.systems, **kwargs)

    def update_plot(self) -> None:
        ket_energy = self.ket.get_energy("GHz")
        energies = [system.get_eigenvalues("GHz") - ket_energy for system in self.systems]
        overlaps = [system.get_eigenbasis().get_overlaps(self.ket) for system in self.systems]

        x_values, xlabel = self.plotwidget._get_x_values_and_label_from_ranges(self.fields)

        self.plotwidget.plot(x_values, energies, overlaps, xlabel=xlabel)

        self.add_short_labels(energies)
        self.plotwidget.add_cursor(x_values[0], energies[0], self.systems[0])

        self.plotwidget.canvas.draw()

    def add_short_labels(
        self,
        energies: Sequence["NDArray[Any]"],
    ) -> None:
        ax = self.plotwidget.canvas.ax
        x_lim = ax.get_xlim()
        ax.set_xlim(x_lim[0] - (x_lim[1] - x_lim[0]) * 0.1, x_lim[1])

        basis0 = self.systems[0].get_eigenbasis()
        corresponding_kets = [basis0.get_corresponding_ket(i) for i in range(basis0.number_of_states)]

        used = set()
        l_dict = {0: "S", 1: "P", 2: "D", 3: "F", 4: "G", 5: "H"}
        for ket, energy in zip(corresponding_kets, energies[0]):
            if "mqdt" not in ket.species:
                short_label = f"{ket.n} {l_dict.get(ket.l, ket.l)}"  # type: ignore
            else:
                short_label = f"{ket.n} L={ket.l:.1f}"
            if short_label in used:
                continue
            used.add(short_label)
            self.plotwidget.canvas.ax.text(x_lim[0], energy, short_label, va="center", ha="right")

    def export_python(self) -> None:
        """Export the current calculation as a Python script."""
        logger.debug("Exporting results as Python script...")

        filename, _ = QFileDialog.getSaveFileName(self, "Save Python Script", "", "Python Files (*.py)")

        if filename:
            filename = filename.removesuffix(".py") + ".py"

            template_path = Path(__file__).parent.parent / "export_templates" / "single_atom.ipynb"
            with open(template_path) as f:
                notebook = nbformat.read(f, as_version=4)

            exporter = PythonExporter(exclude_output_prompt=True, exclude_input_prompt=True)
            content, _ = exporter.from_notebook_node(notebook)

            replacements = self._get_export_replacements()
            for key, value in replacements.items():
                content = content.replace(key, str(value))

            with open(filename, "w") as f:
                f.write(content)

            logger.info(f"Python script saved as {filename}")

    def export_notebook(self) -> None:
        """Export the current calculation as a Jupyter notebook."""
        logger.debug("Exporting results as Jupyter notebook...")

        filename, _ = QFileDialog.getSaveFileName(self, "Save Jupyter Notebook", "", "Jupyter Notebooks (*.ipynb)")

        if filename:
            filename = filename.removesuffix(".ipynb") + ".ipynb"

            template_path = Path(__file__).parent.parent / "export_templates" / "single_atom.ipynb"
            with open(template_path) as f:
                notebook = nbformat.read(f, as_version=4)

            replacements = self._get_export_replacements()
            for cell in notebook.cells:
                if cell.cell_type == "code":
                    source = cell.source
                    for key, value in replacements.items():
                        source = source.replace(key, str(value))
                    cell.source = source

            nbformat.write(notebook, filename)

            logger.info(f"Jupyter notebook saved as {filename}")

    def _get_export_replacements(self) -> dict[str, Any]:
        ket = self.ket_config.get_ket_atom(0)
        qns = self.ket_config.get_quantum_numbers(0)
        restrict_qns = self.basis_config.get_restrict_qns(0)
        fields = self.system_atom_config.get_fields()

        replacements = {
            "SPECIES": f'"{ket.species}"',
            "**QUANTUM_NUMBERS_RESTRICTIONS": ", ".join(f"{k}={v}" for k, v in restrict_qns.items()),
            "**QUANTUM_NUMBERS": ", ".join(f"{k}={v}" for k, v in qns.items()),
            "STEPS": self.system_atom_config.steps_spinbox.value(),
            "EX_MIN": fields["Ex"].min_value,
            "EX_MAX": fields["Ex"].max_value,
            "EY_MIN": fields["Ey"].min_value,
            "EY_MAX": fields["Ey"].max_value,
            "EZ_MIN": fields["Ez"].min_value,
            "EZ_MAX": fields["Ez"].max_value,
            "BX_MIN": fields["Bx"].min_value,
            "BX_MAX": fields["Bx"].max_value,
            "BY_MIN": fields["By"].min_value,
            "BY_MAX": fields["By"].max_value,
            "BZ_MIN": fields["Bz"].min_value,
            "BZ_MAX": fields["Bz"].max_value,
        }

        # Add diagonalization kwargs if needed
        diag_kwargs: dict[str, Any] = {}
        if self.plotwidget.energy_range.isChecked():
            _energies = self.plotwidget.energy_range.values()
            diag_kwargs["energy_range"] = tuple(ket.get_energy("GHz") + v for v in _energies)
            diag_kwargs["energy_unit"] = "GHz"
        if self.plotwidget.fast_mode.isChecked():
            diag_kwargs["diagonalizer"] = "lapacke_evr"
            diag_kwargs["float_type"] = "float32"
        replacements["**DIAGONALIZE_KWARGS"] = ", ".join(f"{k}={v!r}" for k, v in diag_kwargs.items())

        # Handle x-axis values and label
        _, xlabel = self.plotwidget._get_x_values_and_label_from_ranges(fields)
        replacements["X_VALUES"] = xlabel.split(" ")[0]
        replacements["XLABEL"] = f'"{xlabel}"'

        isreal = fields["By"].is_zero() and fields["Ey"].is_zero()
        replacements["PI_DTYPE"] = "real" if isreal else "complex"

        return replacements
