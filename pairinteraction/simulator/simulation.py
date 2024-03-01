"""Defining the simulation class for pairinteraction"""

from typing import Dict, Union

# from pairinteraction.pairinteraction_backend import (
#     BasisAtom,
#     BasisAtomCreator,
#     BasisClassicalLight,
#     ProductBasis,
#     System,
#     SystemAtom,
#     SystemClassicalLight,
#     SystemCombined,
#     SystemWithInteractions,
# )
from pairinteraction.validator.constituents import BaseModelConstituent, ModelAtom, ModelClassicalLight
from pairinteraction.validator.interactions import ModelInteractions
from pairinteraction.validator.model import Model, create_model
from pairinteraction.validator.numerics import ModelNumerics

BasisAtom = (
    BasisAtomCreator
) = (
    BasisClassicalLight
) = (
    ProductBasis
) = (
    System
) = (
    SystemAtom
) = SystemClassicalLight = SystemCombined = SystemWithInteractions = object  # mock to create the documentation


class BaseSimulation:
    """BaseSimulation class for pairinteraction

    This class is used to define a simulation object for pairinteraction.
    It contains the necessary methods to run a simulation and to extract
    the results.
    """

    def __init__(self, settings: Union[Model, Dict]):
        """Initialize the simulation object

        Args:
            settings (dict): Settings for the simulation
        """
        self.settings: Model = create_model(settings)
        self.bases = {}
        self.constituents = {}


class Simulation(BaseSimulation):
    """Simulation class for pairinteraction

    This class is used to define a simulation object for pairinteraction.
    It contains the necessary methods to run a simulation and to extract
    the results.
    """

    def run(self):
        for key, model in self.settings.constituents.items():
            self.finish_model_constituent(model)
            system = self.create_constituent(key, model)
            self.diagonalize_system(system, model)

        model = self.settings.interactions
        self.finish_model_interactions(model)
        system = self.create_system_interactions(model)
        self.diagonalize_system(system, model)

    def finish_model_constituent(self, model: BaseModelConstituent) -> None:
        # TODO set min/max energy if delta energy was given by getting the states of interest energy from the database
        # should all this be done here, or just pass the database, or the states
        # to some function of the constituent model, like:
        # model.set_min_max_energy(database/states)
        pass

    def finish_model_interactions(self, model: ModelInteractions) -> None:
        # TODO set min/max energy if delta energy was given by getting the energies from the diagonalized constituents
        # should this all be done here, or just pass the diagonalized constituents, or the energies
        # to some function of the constituent model, like:
        # model.set_min_max_energy(system_constituent/energies)
        pass

    def create_constituent(self, key: str, model: BaseModelConstituent) -> System:
        """Create the constituent system, i.e. create the basis and the system.

        Returns:
            System: The created system
        """
        if key in ["atom1", "atom2"]:
            basis = self.create_atom_basis(model)
            system = self.create_atom(basis, model)
        elif key in ["classical_light1", "classical_light2"]:
            basis = self.create_classical_light_basis(model)
            system = self.create_classical_light(basis, model)
        self.bases[key] = basis
        self.constituents[key] = system
        return system

    def create_atom_basis(self, model: ModelAtom) -> BasisAtom:
        """Create the basis for an atom.

        This is a simple wrapper around the BasisAtomCreator.
        """
        basis_creator = BasisAtomCreator(model.species)
        basis_creator.restrict_quantum_number_n(model.min_n, model.max_n)
        basis_creator.restrict_quantum_number_nu(model.min_nu, model.max_nu)
        basis_creator.restrict_quantum_number_l(model.min_l, model.max_l)
        basis_creator.restrict_quantum_number_s(model.min_s, model.max_s)
        basis_creator.restrict_quantum_number_j(model.min_j, model.max_j)
        basis_creator.restrict_quantum_number_f(model.min_f, model.max_f)
        basis_creator.restrict_quantum_number_m(model.min_m, model.max_m)
        basis_creator.restrict_energy(model.min_energy, model.max_energy)
        basis_creator.add_additional_states(model.additionally_included_states)
        basis = basis_creator.create()

        # Option 2  # currently not supported, maybe later implemented in python
        # basis_options = {...}
        # basis = StandardBasisAtom(basis_options)  # or **basis_options
        # basis = StandardBasisAtom.from_model(basis_options)  # or **basis_options # .from_options
        return basis

    def create_atom(self, basis: BasisAtom, model: ModelAtom) -> SystemAtom:
        """Create the atom system."""
        # Construct the system
        system = SystemAtom(basis)
        system.set_efield(model.efield)
        system.set_bfield(model.bfield)
        self.set_system_numerics(system, self.settings.numerics)

        # Option 2  # currently not supported, mybe later implemented in python
        # system_options = {...}
        # system = SystemAtom(basis, system_options)  # or analog **system_options or from_model/options
        return system

    def create_classical_light_basis(self, model: ModelClassicalLight) -> BasisClassicalLight:
        raise NotImplementedError

    def create_classical_light(self, basis: BasisClassicalLight, model: ModelClassicalLight) -> SystemClassicalLight:
        raise NotImplementedError

    def create_system_interactions(self, model: ModelInteractions) -> SystemWithInteractions:
        # Construct the basis
        names = sorted(self.constituents.keys())
        bases = [self.bases[k] for k in names]
        basis = self.bases["interactions"] = ProductBasis(*bases, names=names)

        basis.restrict_energy(model.min_energy, model.max_energy)
        basis.set_conserved_total_m(model.conserved_total_m)
        basis.set_conserved_parity_under_inversion(model.conserved_parity_under_inversion)
        basis.set_conserved_parity_under_reflection(model.conserved_parity_under_reflection)
        basis.set_conserved_parity_under_permutation(model.conserved_parity_under_permutation)
        basis.add_additional_states(model.additionally_included_states)  # TODO is this also wanted for combined?

        # TODO is it clever to split SystemCombined creation and adding diagonal energy from the constituents?
        system = SystemCombined(basis)
        for k, constit in self.constituents.items():
            system.add_constituent_energy(k, constit)

        system.set_distance(model.distance)
        system.set_angle(model.angle)

        self.set_system_numerics(system, self.settings.numerics)
        return system

    def set_system_numerics(self, system: System, numerics: ModelNumerics):
        system.set_quantization_axis(numerics.quantization_axis)
        system.set_use_diamagnetism(numerics.use_diamagnetism)
        # TODO etc

    def diagonalize_system(self, system: System, model):
        """Diagonalize the system."""
        # TODO diagonalize options? method, size, min/max energy after, ....
        system.diagonalize()
        # -> this should change the system.basis object to the new eigenbasis
        # is it making a copy of the basis or changing the basis object?

        return system
