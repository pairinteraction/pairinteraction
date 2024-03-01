"""Defining the simulation class for pairinteraction"""

from typing import Dict, List, Union

from pairinteraction.simulator_old.atom import Atom, AtomOne, AtomTwo
from pairinteraction.validator.model import Model, create_model


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

    def create_atom_system(self) -> Atom:
        """Create the atom system

        Returns:
            atom_system: Atom system object
        """
        if self.settings.atom2 is None:
            atom_system = AtomOne(self.settings)
        else:
            atom_system = AtomTwo(self.settings)
        return atom_system

    @staticmethod
    def results_from_atom_system(atom_system: Atom) -> Dict:
        """Extract the results from the atom system

        Args:
            atom_system: Atom system object

        Returns:
            results: Dictionary containing the results of the simulation
        """
        # TODO: save overlaps, ...?
        # TODO: possibility to save the atom_system in results?
        results = {
            "energies": atom_system.energies,
            "settings": atom_system.settings.model_dump(),
        }
        return results

    def run(self):
        raise NotImplementedError("This method should be implemented in the derived class.")


class OneSimulation(BaseSimulation):
    """Simulation class for pairinteraction

    This class is used to define a simulation object for pairinteraction.
    It contains the necessary methods to run a simulation and to extract
    the results.
    """

    def __init__(self, settings: Union[Model, Dict]):
        """Initialize the simulation object

        Args:
            settings (dict): Settings for the simulation
        """
        super().__init__(settings)
        assert (
            self.settings.parameter_range_options.steps == 1
        ), "SimpleSimulation is only for doing one run. Use Simulation for multiple runs."

    def run(self) -> Dict:
        """Run a simple simulation of the settings.

        Returns:
            results: Dictionary containing the results of the simulation
        """
        atom_system = self.create_atom_system()
        results = self.results_from_atom_system(atom_system)
        self.atom_system = atom_system
        self.results = results
        return results


class Simulation(BaseSimulation):
    """Simulation class for pairinteraction

    This class is used to define a simulation object for pairinteraction.
    It contains the necessary methods to run a simulation and to extract
    the results.
    """

    def run(self) -> List[Dict]:
        """Run a simple (not parallelized) simulations of the settings.

        Returns:
            results_list: List of dictionaries containing the results of the simulations.
        """
        results_list = []
        atom_system = self.create_atom_system()
        self.atom_system = atom_system

        for i in range(self.settings.parameter_range_options.steps):
            atom_system.updateParameterStep(i)
            results = self.results_from_atom_system(atom_system)
            results_list.append(results)

        return results_list
