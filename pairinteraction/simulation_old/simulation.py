"""Defining the simulation class for pairinteraction"""

from typing import Dict, List, Union

from pairinteraction.model.model import ModelSimulation
from pairinteraction.simulation_old.atom import Atom, AtomOne, AtomTwo


class BaseSimulation:
    """BaseSimulation class for pairinteraction

    This class is used to define a simulation object for pairinteraction.
    It contains the necessary methods to run a simulation and to extract
    the results.
    """

    def __init__(self, model: Union[ModelSimulation, Dict]):
        """Initialize the simulation object

        Args:
            model (dict): ModelSimulation for the simulation
        """
        self.model: ModelSimulation = ModelSimulation.model_validate(model)

    def create_atom_system(self) -> Atom:
        """Create the atom system

        Returns:
            atom_system: Atom system object
        """
        if self.model.atom2 is None:
            atom_system = AtomOne(self.model)
        else:
            atom_system = AtomTwo(self.model)
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
            "model": atom_system.model.model_dump(),
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

    def __init__(self, model: Union[ModelSimulation, Dict]):
        """Initialize the simulation object

        Args:
            model (dict): ModelSimulation for the simulation
        """
        super().__init__(model)
        assert (
            self.model.parameter_size == 1
        ), "SimpleSimulation is only for doing one run. Use Simulation for multiple runs."

    def run(self) -> Dict:
        """Run a simple simulation of the model.

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
        """Run a simple (not parallelized) simulations of the model.

        Returns:
            results_list: List of dictionaries containing the results of the simulations.
        """
        results_list = []
        atom_system = self.create_atom_system()
        self.atom_system = atom_system

        for i in range(self.model.parameter_size):
            atom_system.updateParameterStep(i)
            results = self.results_from_atom_system(atom_system)
            results_list.append(results)

        return results_list
