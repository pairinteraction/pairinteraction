"""Defining class Atom, AtomOne and AtomTwo.
This classes provide some overhead functionality for calculations with one/two rydberg atoms
and are basically just a wrapper for the pairinteraction.SystemOne and pairinteraction.SystemTwo classes
using the Model class to define the system.
"""
import logging
from pathlib import Path
from typing import List, Tuple, Union

import numpy as np

from pairinteraction import picomplex, pireal
from pairinteraction.simulator_old.exceptions import CppDeleted, CppObjectAlreadyDeleted, QnNotFoundError
from pairinteraction.validator.model import Model
from pairinteraction.validator.states import BaseModelState

logger = logging.getLogger(__name__)


class Atom:
    """This class provides some overhead functionality for calculations with one/two rydberg atoms.

    Parameters
    ----------
    model : dict or Model
        Dictionary (or Model object) containing all relevant parameters for the calculation.

    Attributes
    ----------
    model : Model
        Containing all relevant parameters for the calculation.
    system : pi.(SystemOneReal, SystemOneComplex, SystemTwoReal, SystemTwoComplex)
        This is the actual system.
    cache : pairinteraction.MatrixElementCache
        The MatrixElementCache object used for the calculation.
    basisStates : list of pairinteraction.StateOne or pairinteraction.StateTwo
        List of all basis states in the system, given by system.getMainStates().
    basisQnumbers : list of tuples (n, l, j, m) or ((n1, n2), (l1, l2), (j1, j2), (m1, m2))
        List of all basis qnumbers in the system.
    allStates : list of pairinteraction.StateOne or pairinteraction.StateTwo
        List of all states in the system, given by system.getStates().
    allQnumbers : list of tuples (n, l, j, m) or ((n1, n2), (l1, l2), (j1, j2), (m1, m2))
        List of all qnumbers in the system.
    energies : np.array
        The eigenenergies of the system.
        Calling this property will diagonalize the system if not already done.
    vectors : scipy sparse matrix (len(allStates), len(basisStates))
        The eigenvectors of the system.
        Calling this property will diagonalize the system if not already done.
        Call .toarray() on this to get a dense matrix.
        (or if only a specific row/column is wanted call .getrow(i)/.getcol(i) and then .toarray())
        np.abs(vectors.toarray())**2 corresponds to the overlap of the eigenvectors with all states,
        i.e. the entry [i, j] is the overlap of the eigenvector j with the state i of allStates.
    overlaps : np.array
        Just a shortcut for np.abs(self.vectors.toarray())**2.
    """

    nAtoms = None
    qnumber_basis = ["n", "l", "j", "m"]
    qnumber_types = {"n": int, "l": int, "j": float, "m": float}

    def __init__(self, model: Union[Model, dict]):
        logger.debug("Atom: Initializing Atom object %s.", self.__class__.__name__)
        if not isinstance(model, Model):
            model = Model.model_validate(model)
        self.model = model
        self.s_atom1 = model.atom1
        self.s_atom2 = model.atom2
        self.s_interactions = model.interactions
        self.s_numerics = model.numerics

        self._setdefaultProperties()

    def _setdefaultProperties(self, cache=True, parameter_step=True):
        # properties call by self.property without underscore to ensure, they are created first
        self._system = None
        if cache:
            self._cache = None
        self._basisStates = None
        self._basisQnumbers = None
        self._allStates = None
        self._allQnumbers = None
        self._energies = None
        self._vectors = None
        self._overlaps = None
        if parameter_step:
            self._parameter_step = None

    def __repr__(self):
        return f"{self.__class__.__name__}(...)"

    @property
    def cache(self) -> Union[pireal.MatrixElementCache, picomplex.MatrixElementCache]:
        if self._cache is None:
            path_cache = self.s_numerics.path_cache
            logger.debug("Using cache at %s", path_cache)
            is_real = self.s_atom1.is_real and self.s_atom2.is_real if self.nAtoms == 2 else self.s_atom.is_real
            pirealcomplex = pireal if is_real else picomplex
            if path_cache is None:
                self._cache = pirealcomplex.MatrixElementCache()
            else:
                Path(path_cache).mkdir(exist_ok=True)
                self._cache = pirealcomplex.MatrixElementCache(str(path_cache))

            if self.s_numerics.radial_wavefunction_method == "whittaker":
                self._cache.setMethod(pirealcomplex.WHITTAKER)
            elif self.s_numerics.radial_wavefunction_method == "numerov":
                self._cache.setMethod(pirealcomplex.NUMEROV)
            else:
                raise ValueError(
                    "Unknown MatrixElementCache radial_wavefunction_method "
                    + self.s_numerics.radial_wavefunction_method
                )
            q_defect_db = self.s_numerics.path_quantum_defects_database
            if q_defect_db is not None:
                logger.debug("Using quantum defect database at %s", q_defect_db)
                self._cache.setDefectDB(q_defect_db)
        elif self._cache is CppDeleted:
            raise CppObjectAlreadyDeleted(self.__class__.__name__ + ".cache")
        return self._cache

    @property
    def system(self) -> Union[pireal.SystemOne, picomplex.SystemOne, pireal.SystemTwo, picomplex.SystemTwo]:
        if self._system is None:
            self._createSystem()
        elif self._system is CppDeleted:
            raise CppObjectAlreadyDeleted(self.__class__.__name__ + ".system")
        return self._system

    @property
    def basisStates(self) -> Union[List[pireal.StateOne], List[pireal.StateTwo]]:
        if self._basisStates is None:
            self._createBasis()
        elif self._basisStates is CppDeleted:
            raise CppObjectAlreadyDeleted(self.__class__.__name__ + ".basisStates")
        return self._basisStates

    @property
    def basisQnumbers(self) -> List[tuple]:
        if self._basisQnumbers is None:
            self._createBasis()
        return self._basisQnumbers

    @property
    def allStates(self) -> Union[List[pireal.StateOne], List[pireal.StateTwo]]:
        if self._allStates is None:
            self._createBasis()
        elif self._allStates is CppDeleted:
            raise CppObjectAlreadyDeleted(self.__class__.__name__ + ".allStates")
        return self._allStates

    @property
    def allQnumbers(self) -> List[tuple]:
        if self._allQnumbers is None:
            self._createBasis()
        return self._allQnumbers

    @property
    def energies(self):
        if self._energies is None:
            self._calcEnergies()
        return self._energies

    @property
    def vectors(self):
        if self._vectors is None:
            self._calcEnergies()
        return self._vectors

    @property
    def overlaps(self):
        if self._overlaps is None:
            self._overlaps = np.abs(self.vectors.toarray()) ** 2
        return self._overlaps

    @property
    def parameter_step(self):
        return self._parameter_step

    def _createSystem(self):
        """Creating the actual pi.SystemOne or pi.SystemTwo object.

        Has to be implemented in the subclass, since this really differs for AtomOne or AtomTwo
        """
        if self._system is None:
            raise NotImplementedError
        self.logSystem()

    def _createBasis(self):
        """Creating the basis states and qnumbers and all states and qnumbers.

        allStates is just a list of all states in the system, given by system.getStates().
        basisStates is a list of all basis states in the system, given by system.getMainStates().

        The Qnumbers are just the quantum numbers of the states, for a single state,
        there are given dependent wether we have nAtoms == 1 or 2:
        StateOne Qnumbers: (n, l, j, m)
        StateTwo Qnumbers: ((n1, n2), (l1, l2), (j1, j2), (m1, m2))
        """
        try:
            self._allStates = list(self.system.getStates())
            self._basisStates = list(self.system.getMainStates())
        except RuntimeError as e:
            logger.warning("%s WARNING: %s", self, e)
            self._allStates = [] if self._allStates is None else self._allStates
            self._basisStates = [] if self._basisStates is None else self._basisStates

        self._basisQnumbers = [self.stateToQn(s) for s in self._basisStates]
        self._allQnumbers = [self.stateToQn(s) for s in self._allStates]

    def stateToQn(self, state):
        """Get the quantum numbers from a state object."""
        qn = (state.getN(), state.getL(), state.getJ(), state.getM())
        return self.formatQnumber(qn)

    def _calcEnergies(self):
        """Calculating all energies for a given system."""
        system = self.system
        logger.debug("Calculating energies for %s", type(system).__name__)
        system.diagonalize(self.s_numerics.diagonalize_threshold)
        self._energies = np.real_if_close(system.getHamiltonian().diagonal())  # np.array (len(basisStates))
        self._vectors = system.getBasisvectors()  # sparse matrix (len(allStates), len(basisStates))
        self._createBasis()  # just in case something reordered # TODO should not happen, but does it happen?

    def logSystem(self, loglevel="debug"):
        log = {
            "debug": logger.debug,
            "info": logger.info,
            "warning": logger.warning,
            "print": print,
        }[loglevel]
        log(type(self).__name__)
        log("Number of basis vectors: %d", len(self.basisStates))
        log("Number of involved states: %d", len(self.allStates))

    def formatQnumber(self, qn) -> tuple:
        raise NotImplementedError("This has to be implemented in the subclass")

    def getQnIndex(self, qn) -> int:
        """Finding the index of a given quantum number in the list of basisQnumbers.

        Parameters
        ----------
        qn : qn-like tuple or list (see self.formatQnumber)
            Quantum numbers (e.g. (n, l, j, m) for AtomOne or ((n1, n2), (l1, l2), (j1, j2), (m1, m2)) for two atoms.

        Returns
        -------
        index : int
            Index of the given quantum number in the list of basis states.
        """
        qn = self.formatQnumber(qn)
        try:
            return self.basisQnumbers.index(qn)
        except ValueError as exc:
            if qn in self.allQnumbers:
                msg = (
                    f"getQnIndex: {qn} is not in basisQnumbers, but in allQnumbers,"
                    "if you wanted to get the index of qn in allQnumbers you should call getQnIndexAll."
                )
                raise QnNotFoundError(qn, useAll=False, msg=msg) from exc
            raise QnNotFoundError(qn, useAll=False) from exc

    def getQnIndexAll(self, qn) -> int:
        """Finding the index of a given quantum number in the list of allQnumbers.

        Parameters
        ----------
        qn : qn-like tuple or list (see self.formatQnumber)
            Quantum numbers (e.g. (n, l, j, m) for AtomOne or ((n1, n2), (l1, l2), (j1, j2), (m1, m2)) for two atoms.

        Returns
        -------
        index : int
            Index of the given quantum number in the list of basis states.
        """
        qn = self.formatQnumber(qn)
        try:
            return self.allQnumbers.index(qn)
        except ValueError as exc:
            raise QnNotFoundError(qn, useAll=True) from exc

    def getQnState(self, qn) -> Union[pireal.StateOne, pireal.StateTwo]:
        """Return the state corresponding to the given quantum number.

        Parameters
        ----------
        qn : qn-like tuple or list (see self.formatQnumber)
            Quantum numbers (e.g. (n, l, j, m) for AtomOne or ((n1, n2), (l1, l2), (j1, j2), (m1, m2)) for two atoms.

        Returns
        -------
        state : pairinteraction.StateOne or pairinteraction.StateTwo
            The state corresponding to the given quantum number.
        """
        # TODO if not in allStats should we just newly construct the state?
        return self.allStates[self.getQnIndexAll(qn)]

    def getStateEnergy(self, state: Union[pireal.StateOne, BaseModelState, pireal.StateTwo], **kwargs) -> float:
        """Return the eigenenergy of the state closest to the given state.

        Parameters
        ----------
        state : pairinteraction.StateOne, pairinteraction.StateTwo or instance of BaseModelState
            The state to get the energy from.

        Returns
        -------
        energy : float
            The energy of the state.
        """
        if isinstance(state, BaseModelState):
            state = state.state
        qn = self.stateToQn(state)
        return self.getQnEnergy(qn, **kwargs)

    def getQnEnergy(self, qn, *, ignoreAmbigous=False) -> float:
        """Return the energy corresponding to the given quantum number.
        Corresponding here means the energies which eigenvectors have the largest overlap with the given qn.

        Parameters
        ----------
        qn : qn-like tuple or list (see self.formatQnumber)
            Quantum numbers (e.g. (n, l, j, m) for AtomOne or ((n1, n2), (l1, l2), (j1, j2), (m1, m2)) for two atoms.
        ignoreAmbigous : bool, optional
            If the energy is ambigous still return the most likely energy instead of np.nan.
            Default is False.

        Returns
        -------
        energy : float
            The energy corresponding to the given states.
            If it is not possible to uniquely determine the energy, np.nan is returned.
        """
        try:
            qnOverlaps = self.getQnOverlaps(qn)
        except QnNotFoundError:
            logger.warning("%s.getQnEnergy: qn %s not found, returning nan as energy.", self, qn)
            return np.nan

        inds = np.argsort(qnOverlaps)[::-1]
        sortedEnergies = self.energies[inds]
        sortedOverlaps = qnOverlaps[inds]

        if sortedOverlaps[0] > 0.55 or ignoreAmbigous:
            return sortedEnergies[0]

        sameEnergies = np.isclose(sortedEnergies, sortedEnergies[0], rtol=1e-05, atol=1e-08)
        if sum(sameEnergies) > 1:
            summedOverlaps = np.sum(sortedOverlaps[sameEnergies])
            if summedOverlaps > 0.55:
                return sortedEnergies[0]

        logger.warning(
            "%s.getQnEnergy: overlap of %s is ambigous: overlaps %s with energies %s",
            self,
            qn,
            sortedOverlaps[:5],
            sortedEnergies[:5],
        )
        return np.nan

    def getQnsEnergies(self, qns) -> np.array:
        """Return the energies corresponding to the given quantum numbers.
        See getQnEnergy for more details.
        """
        return np.array([self.getQnEnergy(qn) for qn in qns])

    def getQnOverlaps(self, qn) -> np.array:
        """Return the overlaps of the eigenstates with the state specified by the given quantum number.

        Parameters
        ----------
        qn : qn-like tuple or list (see self.formatQnumber)
            Quantum numbers (e.g. (n, l, j, m) for AtomOne or ((n1, n2), (l1, l2), (j1, j2), (m1, m2)) for two atoms.

        Returns
        -------
        overlaps : np.array
            The overlaps of the eigenstates with the specified state.
        """
        ind = self.getQnIndexAll(qn)
        try:
            return np.abs(self.vectors.getrow(ind).toarray().flatten()) ** 2
        except QnNotFoundError:
            return np.zeros(len(self.basisStates))

    def updateParameterStep(self, step: int) -> bool:
        if self._parameter_step == step:
            return False
        parameter_range_options = self.model.parameter_range_options
        if not 0 <= step < parameter_range_options.steps:
            raise ValueError(f"Invalid step {step} for AtomOne.updateParameterStep")
        self._parameter_step = step
        return True

    def copySystem(self):
        """Returning a copy of the system."""
        return type(self.system)(self.system)

    def deleteCppObjects(self):
        """Delete all cpp objects, usually you want to call this, before saving the results
        to a pkl file to make the size much smaller and prevent weird bugs due to the cpp objects.
        """
        # Never delete _cache without also deleting the _system!
        # Also first delete _system, then _cache
        # If deleting _system and _cache first copy the _vectors object (scipy.sparse._csc.csc_matrix),
        # otherwise it will be corrupt after deleting the _system/_cache
        if self._vectors is not None:
            self._vectors = self._vectors.copy()
        for k in [
            "_system",
            "_cache",
            "_basisStates",
            "_allStates",
        ]:
            if hasattr(self, k):
                delattr(self, k)
                setattr(self, k, CppDeleted)


class AtomOne(Atom):
    nAtoms = 1

    def __init__(self, model: Union[Model, dict], iAtom: int = 1):
        super().__init__(model)

        self.iAtom = int(iAtom)
        self.s_atom = getattr(self.model, f"atom{self.iAtom}")

    def _createSystem(self):
        """Creating the actual pi.SystemOne."""
        s_atom = self.s_atom
        SystemOne = pireal.SystemOne if s_atom.is_real else picomplex.SystemOne
        self._system = system = SystemOne(s_atom.species, self.cache)

        for q in ["n", "l", "j", "m"]:
            minQ, maxQ = getattr(s_atom, "min_" + q), getattr(s_atom, "max_" + q)
            if minQ is not None and maxQ is not None:
                restrict = getattr(system, "restrict" + q.capitalize())
                restrict(minQ, maxQ)
        minEnergy, maxEnergy = s_atom.min_energy, s_atom.max_energy
        if minEnergy is not None and maxEnergy is not None:
            system.restrictEnergy(minEnergy, maxEnergy)

        system.enableDiamagnetism(self.s_numerics.use_diamagnetism)
        self.setBfield()
        self.setEfield()

        super()._createSystem()

    def setBfield(self):
        """Setting the magnetic field for the system."""
        bfield = [getattr(self.s_atom, f"bfield_{c}").get_value(self.parameter_step) for c in "xyz"]
        if bfield != [0, 0, 0]:
            self.system.setBfield(bfield)

    def setEfield(self):
        """Setting the electric field for the system."""
        efield = [getattr(self.s_atom, f"efield_{c}").get_value(self.parameter_step) for c in "xyz"]
        if efield != [0, 0, 0]:
            self.system.setEfield(efield)

    @classmethod
    def formatQnumber(cls, qn) -> Tuple[int, int, float, float]:
        """Format the qn to tuples like (n, l, j, m).

        Parameters
        ----------
        qn : tuple or list
            Like (n, l, j, m)

        Returns
        -------
        qn : tuple
            Tuple (n, l, j, m)

        """
        shape = np.shape(qn)
        if shape == (4,):
            pass
        else:
            if len(shape) > 1:
                raise ValueError(
                    "Can't recognize shape of qn. If you want to format a list of qn use formatQnumberList instead."
                )
            raise ValueError("Can't recognize shape of qn.")

        qn = tuple(cls.qnumber_types[q](qq) for q, qq in zip(cls.qnumber_basis, qn))
        return qn

    def updateParameterStep(self, step: int) -> bool:
        """Updating the paramter step and set all attributes, that have to be computed newly to None.

        Returns
        -------
        updated : bool
            Whether the system was updated.
        """
        updated = super().updateParameterStep(step)
        if not updated:
            return False
        parameters = self.model.parameter_range_options.parameters

        if any(k in parameters for k in ["efield_x", "efield_y", "efield_z", "bfield_x", "bfield_y", "bfield_z"]):
            self._energies, self._vectors, self._overlaps = None, None, None
            self.setEfield()
            self.setBfield()
            return True
        return False


class AtomTwo(Atom):
    nAtoms = 2

    def __init__(self, model: Union[Model, dict]):
        super().__init__(model)

        # properties call by self.property without underscore to ensure, they are created first
        self._atom1 = None
        self._atom2 = None

    @property
    def atom1(self) -> AtomOne:
        return self.getAtom(1)

    @property
    def atom2(self) -> AtomOne:
        return self.getAtom(2)

    @property
    def atoms(self) -> List[AtomOne]:
        if self.useSameAtom:
            return [self.atom1]
        return [self.atom1, self.atom2]

    @property
    def useSameAtom(self) -> bool:
        return isinstance(self.s_atom2, str) and self.s_atom2 == "atom1"

    def getAtom(self, iAtom) -> AtomOne:
        """Getting the AtomOne object for atom 1 or 2."""
        if iAtom == 2 and self.useSameAtom:
            return self.getAtom(1)

        if getattr(self, f"_atom{iAtom}") is None:
            atom = AtomOne(self.model, iAtom=iAtom)
            setattr(self, f"_atom{iAtom}", atom)
        return getattr(self, f"_atom{iAtom}")

    def _createSystem(self):
        """Creating the actual pi.SystemTwo, based on the AtomOne.system's."""
        for atom in self.atoms:
            atom.energies  # ensure the system is diagonalized first  # noqa: B018

        is_real = self.s_atom1.is_real and self.s_atom2.is_real
        SystemTwo = pireal.SystemTwo if is_real else picomplex.SystemTwo
        self._system = system = SystemTwo(self.atom1.system, self.atom2.system, self.cache)

        s_interactions = self.s_interactions
        s_numerics = self.s_numerics

        if s_interactions.use_delta_energy_after_fields:
            energies = []
            for csoi in s_interactions.combined_states_of_interest:
                energy = 0
                for k, s in csoi.items():
                    constituent = getattr(self, k, None)
                    if constituent is None:
                        raise ValueError(f"Constituent {k} not found in AtomTwo.")
                    energy += constituent.getQnEnergy((s.n, s.l, s.j, s.m), ignoreAmbigous=True)
                energies.append(energy)
            delta = s_interactions.delta_energy
            min_energy, max_energy = min(energies) - delta, max(energies) + delta
        else:
            min_energy, max_energy = s_interactions.min_energy, s_interactions.max_energy

        if min_energy is not None and max_energy is not None:
            if any(np.isnan([min_energy, max_energy])):
                raise ValueError("AtomTwo min/max_energy has nan values.")
            system.restrictEnergy(min_energy, max_energy)

        if s_numerics.precision_sparse_matrices is not None:
            system.setMinimalNorm(s_numerics.precision_sparse_matrices)
        if s_interactions.conserved_total_m is not None:
            system.setConservedMomentaUnderRotation(s_interactions.conserved_total_m)
        if s_numerics.rydberg_rydberg_multipole_order is not None:
            system.setOrder(s_numerics.rydberg_rydberg_multipole_order)

        self.setSymmetry("all")
        self.setDistance()
        self.setAngle()

        super()._createSystem()

    def setSymmetry(self, sym: str):
        """Setting the conserved parity under a certain symmetry."""
        if sym == "all":
            for s in ["inversion", "permutation", "reflection"]:
                self.setSymmetry(s)
            return

        value = getattr(self.s_interactions, f"conserved_parity_under_{sym}").get_value(self.parameter_step)
        if value is None:
            return
        elif value not in [-1, 1]:
            raise ValueError(f"Unknown value {value} for conserved_parity_under_{sym}.")
        logger.debug("Setting %s to %s", sym, value)
        oddeven = pireal.ODD if value == -1 else pireal.EVEN
        getattr(self.system, "setConservedParityUnder" + sym.capitalize())(oddeven)

    def setDistance(self):
        """Setting the distance for the system."""
        distance = self.s_interactions.distance.get_value(self.parameter_step)
        if distance is not None and not np.isinf(distance):
            self.system.setDistance(distance)

    def setAngle(self):
        """Setting the angle for the system."""
        angle = self.s_interactions.angle.get_value(self.parameter_step)
        if angle is not None:
            self.system.setAngle(angle)

    @classmethod
    def formatQnumber(cls, qn) -> Tuple[Tuple[int, int], Tuple[int, int], Tuple[float, float], Tuple[float, float]]:
        """Format the qn to tuples like ((n1, n2), (l1, l2), (j1, j2), (m1, m2)).

        Parameters
        ----------
        qn : tuple or list
            Takes care of the different possible formats:
            - ((n1, l1, j1, m1), (n2, l2, j2, m2))
            - (n1, l1, j1, m1, n2, l2, j2, m2)
            - ((n1, n2), (l1, l2), (j1, j2), (m1, m2)) (this is the recommended format)

        Returns
        -------
        qn : tuple
            Tuple like: ((n1, n2), (l1, l2), (j1, j2), (m1, m2))

        """
        shape = np.shape(qn)
        if shape == (8,):
            qn = tuple((qn[i], qn[i + 4]) for i in range(4))
        elif shape == (4, 2):
            pass
        elif shape == (2, 4):
            qn = tuple(tuple(q) for q in np.transpose(qn))
        else:
            if len(shape) > 2:
                raise ValueError(
                    "Can't recognize shape of qn. If you want to format a list of qn use formatQnumberList instead."
                )
            raise ValueError("Can't recognize shape of qn.")
        qn = tuple(tuple(cls.qnumber_types[q](qq) for qq in qqs) for q, qqs in zip(cls.qnumber_basis, qn))
        return qn

    def updateParameterStep(self, step: int) -> bool:
        """Updating the paramter step and set all attributes, that have to be computed newly to None.

        Returns
        -------
        updated : bool
            Whether the system was updated.
        """
        updated = super().updateParameterStep(step)
        if not updated:
            return False

        atoms_updated = False
        for atom in self.atoms:
            atoms_updated = atom.updateParameterStep(step) or atoms_updated

        parameters = self.model.parameter_range_options.parameters
        if atoms_updated or any(
            k in parameters for k in [f"conserved_parity_under_{x}" for x in ["inversion", "permutation", "reflection"]]
        ):
            self._setdefaultProperties(cache=False, parameter_step=False)
            self._createSystem()
            self._createBasis()
            return True
        elif any(k in parameters for k in ["distance", "angle"]):
            self._energies, self._vectors, self._overlaps = None, None, None
            self.setDistance()
            self.setAngle()
            return True
        # TODO onlySameTrafo

        return False

    def deleteCppObjects(self):
        """Delete all cpp objects, usually you want to call this, before saving the results
        to a pkl file to make the size much smaller and prevent weird bugs due to the cpp objects.
        """
        for k in ["_atom1", "_atom2"]:
            if getattr(self, k, None) is not None:
                getattr(self, k).deleteCppObjects()
        super().deleteCppObjects()
