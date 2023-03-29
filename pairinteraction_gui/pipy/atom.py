"""Defining class Atom, AtomOne and AtomTwo.
This classes provide some overhead functionality for calculations with one/two rydberg atoms
and are basically just a wrapper for the pairinteraction.SystemOne and pairinteraction.SystemTwo classes
using the Config class to define the system.
"""
import logging
import os

import numpy as np
from pipy import Config
from pipy import picomplex
from pipy import pireal

logger = logging.getLogger(__name__)


class Atom:
    """This class provides some overhead functionality for calculations with one/two rydberg atoms.

    Parameters
    ----------
    config : dict or Config
        Dictionary (or Config object) containing all relevant parameters for the calculation.

    Attributes
    ----------
    config : Config
        Containing all relevant parameters for the calculation.
    system : pairinteraction.SystemOne or pairinteraction.SystemTwo
        This is the actual system, either a pairinteraction.SystemOne or pairinteraction.SystemTwo object,
        depending on the subclass.
    basisStates : list of pairinteraction.StateOne or pairinteraction.StateTwo
        List of all basis states in the system, given by system.getMainStates().
        The states are given as pairinteraction.StateOne or pairinteraction.StateTwo objects.
    basisQunumbers : list of tuples (N, L, J, M)
        List of all basis states in the system.
        The states are given as tuples of the quantum numbers N, L, J, M.
    """

    nAtoms = None
    Qunumber_basis = ["N", "L", "J", "M"]

    def __init__(self, config):
        if isinstance(config, dict):
            config = Config(config)
        self.config = config

        # properties call by self.property without underscore to ensure, they are created first
        self._system = None
        self._basisStates = None
        self._basisQunumbers = None
        self._allStates = None
        self._allQunumbers = None

        # this will be set when doing calcEnergies and reset to None if the config changes
        self.energies = None
        self.vectors = None
        self.overlaps = None

    def __repr__(self):
        config = self.config
        return f"{self.__class__.__name__}({config.species()}, ...)"

    @property
    def system(self):
        if self._system is None:
            self._createSystem()
        return self._system

    @property
    def basisStates(self):
        if self._basisStates is None:
            self._createBasis()
        return self._basisStates

    @property
    def basisQunumbers(self):
        if self._basisQunumbers is None:
            self._createBasis()
        return self._basisQunumbers

    @property
    def allStates(self):
        if self._allStates is None:
            self._createBasis()
        return self._allStates

    @property
    def allQunumbers(self):
        if self._allQunumbers is None:
            self._createBasis()
        return self._allQunumbers

    def _createBasis(self):
        self._basisStates = list(self.system.getMainStates())
        self._allStates = list(self.system.getStates())
        basisQunumbers = [(s.getN(), s.getL(), s.getJ(), s.getM()) for s in self._basisStates]
        allQunumbers = [(s.getN(), s.getL(), s.getJ(), s.getM()) for s in self._allStates]
        if self.nAtoms == 2:
            basisQunumbers = [tuple(tuple(q) for q in qn) for qn in basisQunumbers]
            allQunumbers = [tuple(tuple(q) for q in qn) for qn in allQunumbers]
        self._basisQunumbers = basisQunumbers
        self._allQunumbers = allQunumbers

    def copy(self):
        raise NotImplementedError("This has to be implemented in the subclass")

    def findQunumberIndex(self, qunumber):
        """Finding the index of a given quantum number in the list of basisStates.

        Parameters
        ----------
        qunumber : tuple
            Tuple of quantum numbers (N, L, J, M) or ((N1, L1, J1, M1), (N2, L2, J2, M2)) for two atoms.

        Returns
        -------
        index : int
            Index of the given quantum number in the list of basis states.
        """
        if len(qunumber) == len(self.Qunumber_basis):
            pass
        elif len(qunumber) == self.nAtoms and len(qunumber[0]) == len(self.Qunumber_basis):
            logger.warning(
                f"qunumber {qunumber} is given as (N1, L1, J1, M1), (N2, L2, J2, M2), "
                "but should be ((N1, N2), (L1, L2), (J1, J2), (M1, M2), converting it now."
            )
            qunumber = tuple(zip(*qunumber))
        else:
            raise ValueError(f"qunumber {qunumber} has wrong length")

        qunumber = tuple(qunumber)
        if self.nAtoms > 1:
            qunumber = tuple(tuple(q) for q in qunumber)
        return self.basisQunumbers.index(qunumber)

    def _createSystem(self):
        """Creating the actual pi.SystemOne or pi.SystemTwo object.

        Has to be implemented in the subclass, since this really differs for AtomOne or AtomTwo
        """
        if self._system is None:
            raise NotImplementedError
        self.logSystem()

    def logSystem(self):
        logger.debug(type(self).__name__)
        system = self.system

        logger.debug("Number of basis vectors: %s, %s", system.getNumBasisvectors(), len(system.getMainStates()))
        # logger.debug("Basisvectors = MainStates keys [:5]: %s", [s.getKey() for s in system.getMainStates()][:5])

        logger.debug("Number of involved states: %s, %s", system.getNumStates(), len(system.getStates()))
        # logger.debug("States keys [:5]: %s", [s.getKey() for s in system.getStates()][:5])

    def copySystem(self):
        """Returning a copy of the system."""
        return self.system.__class__(self.system)

    def calcEnergies(self, sort=False, calc_overlaps=False):
        """Calculating all energies for a given system.

        Returns
        -------
        energies : np.array
            Sorted array of energies
        overlaps : np.array
            Array of overlaps of allStates with the eigenstates of the corresponding energies.
            The shape is (len(self.allStates), len(energies))
        """
        if self.energies is not None:
            logger.warning("calcEnergies was already called, sure you want to recalculate the energies?")
        system = self.system
        logger.debug("Calculating energies for %s", type(system).__name__)
        system.diagonalize(self.config.diagonalizeThreshold())
        self.energies = np.real_if_close(system.getHamiltonian().diagonal())
        self.vectors = system.getBasisvectors()  # sparse matrix (getNumStates, getNumBasisvectors)
        self._createBasis()  # just in case something reordered

        if sort:
            self.sorted_inds = np.argsort(self.energies)
            self.energies = self.energies[self.sorted_inds]
            self.vectors = self.vectors[:, self.sorted_inds]
        else:
            self.sorted_inds = None

        if calc_overlaps:  # this is slow, so only do it if necessary
            self.overlaps = np.array([system.getOverlap(state) for state in self.allStates])
            if sort:
                self.overlaps = self.overlaps[:, self.sorted_inds]

    def getEnergiesStates(self, qunumbers):
        """Just a wrapper around calc_energies, to give back only the energies corresponding to the given states.
        (Corresponding here means the energies which eigenvectors have the largest overlap with the given states)

        Parameters
        ----------
        qunumbers : list of qunumber
            The quantum numbers of the states the energy should be calculated for.
            The tuple should be of the form (n, l, j, m) for AtomOne
            and (n1, l1, j1, m1, n2, l2, j2, m2) for AtomTwo.

        Returns
        -------
        energies : np.array
            The energies corresponding to the given states.
        """
        if self.energies is None:
            raise ValueError("Energies have not been calculated yet. Call calcEnergies() first.")

        if self.overlaps is not None:
            overlaps_all = self.overlaps
            inds = [self.findQunumberIndex(qunumber) for qunumber in qunumbers]
            overlaps = overlaps_all[inds, :]
        elif self.nAtoms == 1:
            pi = pireal if self.config.isReal() else picomplex
            qunumbers = [(int(q[0]), int(q[1]), float(q[2]), float(q[3])) for q in qunumbers]
            states = [pi.StateOne(self.config.species(), *q) for q in qunumbers]
            overlaps = np.array([self.system.getOverlap(state) for state in states])
            if self.sorted_inds is not None:
                overlaps = overlaps[:, self.sorted_inds]
        elif self.nAtoms == 2:
            pi = pireal if self.config.isReal() else picomplex
            qunumbers = [
                ([int(q[0]), int(q[4])], [int(q[1]), int(q[5])], [float(q[2]), float(q[6])], [float(q[3]), float(q[7])])
                for q in qunumbers
            ]
            species = [self.config.species(1), self.config.species(2)]
            states = [pi.StateTwo(species, *q) for q in qunumbers]
            overlaps = np.array([self.system.getOverlap(state) for state in states])
            if self.sorted_inds is not None:
                overlaps = overlaps[:, self.sorted_inds]

        energies = []
        for i, overlap in enumerate(overlaps):
            inds = np.argsort(overlap)[::-1]
            energies.append(self.energies[inds[0]])

            overlap = overlap[inds]
            if overlap[0] > 0.6 or overlap[0] / overlap[1] > 1.5:
                pass
            elif sum(overlap[:2]) > 2.1 / 3:
                if abs(self.energies[inds[0]] - self.energies[inds[1]]) > 1e-4:
                    logger.warning(
                        "getEnergiesStates: overlap of %s is ambigous and energies differ by %s",
                        qunumbers[i],
                        self.energies[inds[0]] - self.energies[inds[1]],
                    )
            else:
                logger.warning(
                    "getEnergiesStates: overlap of %s is ambigous: overlaps %s with energies %s",
                    qunumbers[i],
                    overlap[:5],
                    self.energies[inds[:5]],
                )
        return np.array(energies)

    def getEnergyState(self, qunumber):
        """Just a wrapper around calc_states_energies for calculating the energy of a single qunumber."""
        return self.getEnergiesStates([qunumber])[0]

    def updateFromParams(self, params):
        raise NotImplementedError("This has to be implemented in the subclass")

    def getCache(self):
        if not hasattr(self, "_cache") or self._cache is not None:
            pathCache = self.config.pathCache()
            os.makedirs(pathCache, exist_ok=True)
            logger.debug("Using cache at %s", pathCache)
            pi = pireal if self.config.isReal() else picomplex
            if self.config.method() == "WHITTAKER":
                pi.MatrixElementCache.setMethod(pi.WHITTAKER)
            self._cache = pi.MatrixElementCache(pathCache)
        return self._cache


class AtomOne(Atom):
    nAtoms = 1

    @property
    def system(self) -> pireal.SystemOne:
        return super().system

    @property
    def basisStates(self) -> "list[pireal.StateOne]":
        return super().basisStates

    def copy(self):
        """Returning a copy of the AtomOne without any energies, overlaps or basisStates.
        Basically a copy of the config and a copy of the system.
        """
        config = self.config.deepCopy()
        new = AtomOne(config)
        new._system = self.copySystem()
        return new

    def _createSystem(self):
        """Creating the actual pi.SystemOne."""
        config = self.config
        pi = pireal if config.isReal() else picomplex
        system = pi.SystemOne(config.species(), self.getCache())

        for Q in ["N", "L", "J", "M"]:
            minQ, maxQ = config.restrictQnSingle(Q)
            if minQ is not None and maxQ is not None:
                restrict = getattr(system, "restrict" + Q)
                restrict(minQ, maxQ)
        minEnergy, maxEnergy = config.restrictEnergySingle(pi=pi)
        if minEnergy is not None and maxEnergy is not None:
            system.restrictEnergy(minEnergy, maxEnergy)

        if config.conserveMomentaSingle():
            system.setConservedMomentaUnderRotation(config.momentaSingle())
        if config.Efield() != [0, 0, 0]:
            system.setEfield(config.Efield())
        if config.Bfield() != [0, 0, 0]:
            system.setBfield(config.Bfield())
        system.enableDiamagnetism(config.diamagnetism())

        self._system = system
        super()._createSystem()

    def updateFromParams(self, new):
        """Updating the system from given parameters.

        Parameters
        ----------
        params : dict
            The parameters to update the system with.
            Should only contain the following keys:
            Ex, Ey, Ez, Bx, By, Bz

        Returns
        -------
        updated : bool
            Whether the system was updated.
        """
        if any(k not in ["Ex", "Ey", "Ez", "Bx", "By", "Bz"] for k in new):
            raise ValueError("Invalid key in AtomOne updateFromParams")

        config = self.config
        Efield = [new.get("Ex", config.Ex()), new.get("Ey", config.Ey()), new.get("Ez", config.Ez())]
        Bfield = [new.get("Bx", config.Bx()), new.get("By", config.By()), new.get("Bz", config.Bz())]
        update = {
            "Efield": Efield != config.Efield(),
            "Bfield": Bfield != config.Bfield(),
        }

        if update["Efield"]:
            config.setEfield(Efield, "")
            self.system.setEfield(config.Efield())
        if update["Bfield"]:
            config.setBfield(Bfield, "")
            self.system.setBfield(config.Bfield())

        if any(update.values()):
            self.energies, self.overlaps, self.vectors = None, None, None
            return True
        return False


class AtomTwo(Atom):
    nAtoms = 2

    def __init__(self, config):
        super().__init__(config)

        # properties call by self.property without underscore to ensure, they are created first
        self._atom1 = None
        self._atom2 = None

    @property
    def system(self) -> pireal.SystemTwo:
        return super().system

    @property
    def basisStates(self) -> "list[pireal.StateTwo]":
        return super().basisStates

    @property
    def atom1(self) -> AtomOne:
        if self._atom1 is None:
            self._createAtom(1)
        return self._atom1

    @property
    def atom2(self) -> AtomOne:
        if self.config.useSameAtom():
            return self.atom1
        if self._atom2 is None:
            self._createAtom(2)
        return self._atom2

    @property
    def atoms(self) -> "list[AtomOne]":
        if self.config.useSameAtom():
            return [self.atom1]
        return [self.atom1, self.atom2]

    def copy(self):
        """Returning a copy of the AtomTwo without any energies, overlaps or basisStates.
        Basically a copy of the config, the system and the atoms.
        """
        config = self.config.deepCopy()
        new = AtomTwo(config)
        if self._atom1 is not None:
            new._atom1 = self.atom1.copy()
        if self._atom2 is not None and not self.config.useSameAtom():
            new._atom2 = self.atom2.copy()
        new._system = self.copySystem()
        return new

    def _createAtom(self, iAtom):
        """Creating the AtomOne objects for atom 1 or 2."""
        config = self.config.shallowCopy()
        config.force_nAtoms = 1
        config.iAtom = iAtom
        atom = AtomOne(config)
        setattr(self, f"_atom{iAtom}", atom)

    def _createSystem(self):
        """Creating the actual pi.SystemTwo, based on the AtomOne.system's."""
        config = self.config
        pi = pireal if config.isReal() else picomplex
        for atom in self.atoms:
            atom.calcEnergies()
        system = pi.SystemTwo(self.atom1.system, self.atom2.system, self.getCache())

        for q in ["N", "L", "J", "M"]:
            minq, maxq = config.restrictQnPair(q)
            if minq is not None and maxq is not None:
                raise NotImplementedError("Restricting QN for SystemTwo is not implemented yet.")
        minEnergy, maxEnergy = config.restrictEnergyPair(atom1=self.atom1, atom2=self.atom2)
        if minEnergy is not None and maxEnergy is not None:
            system.restrictEnergy(minEnergy, maxEnergy)

        if config.minimalNorm() is not None:
            system.setMinimalNorm(config.minimalNorm())
        if config.conserveMomentaPair():
            system.setConservedMomentaUnderRotation(config.momentaPair())
        if config.angle() is not None:
            system.setAngle(config.angle())
        if not np.isinf(config.distance()) and config.distance() is not None:
            system.setDistance(config.distance())
        if config.order() is not None:
            system.setOrder(config.order())

        for sym in ["inversion", "permutation", "reflection"]:
            value = config.symmetry(sym)
            if value is None:
                continue
            logger.debug("Setting %s to %s", sym, value)
            oddeven = pi.ODD if value == "ODD" else pi.EVEN
            getattr(system, "setConservedParityUnder" + sym.capitalize())(oddeven)

        self._system = system
        super()._createSystem()

    def updateFromParams(self, _new):
        """Updating the system from given parameters.

        Parameters
        ----------
        params : dict
            The parameters to update the system with.
            Should only contain the following keys:
                For AtomOne: Ex, Ey, Ez, Bx, By, Bz
                For AtomTwo: distance, angle, inversion, permutation, reflection

        Returns
        -------
        updated : bool
            Whether the system was updated.
        """
        # TODO onlySameTrafo
        new = _new.copy()

        atom_new = {}
        for k in list(new.keys()):
            if k.startswith("atom."):
                atom_new[k[5:]] = new.pop(k)
            elif k in ["Ex", "Ey", "Ez", "Bx", "By", "Bz"]:
                atom_new[k] = new.pop(k)

        if any(k not in ["distance", "angle", "inversion", "permutation", "reflection"] for k in new):
            raise ValueError(f"Invalid key {k} in AtomTwo updateFromParams")

        config = self.config
        update = {
            **{f"atom{i}": atom.updateFromParams(atom_new) for i, atom in enumerate(self.atoms)},
        }
        for k in ["angle", "distance", "inversion", "permutation", "reflection"]:
            if k in new:
                update[k] = new[k] != getattr(config, k)()

        if update.get("angle", False):
            config.setAngle(new["angle"])
            self.system.setAngle(config.angle())
        if update.get("distance", False):
            config.setDistance(new["distance"])
            self.system.setDistance(config.distance())
        for sym in ["inversion", "permutation", "reflection"]:
            if update.get(sym, False):
                config.setSymmetry(sym, new[sym])

        if any(update.get(k, False) for k in ["atom1", "atom2"]):
            # for now we have to build the system again, because internal SystemTwo makes a copy of the SystemOne's ...
            # and therefore by updatind the SystemOne's it does not update the SystemTwo
            logger.warning(
                "AtomTwo updating the E/B field. This will update the AtomOne's, "
                "and has to recalculate the system (but will not change the basis!). Sure you want to do this?"
            )
        if any(update.get(k, False) for k in ["inversion", "permutation", "reflection"]):
            logger.warning(
                "AtomTwo updating a symmetry. This will update the basis and recalculates the system! "
                "Sure you want to do this?"
            )

        if any(update.get(k, False) for k in ["atom1", "atom2", "inversion", "permutation", "reflection"]):
            self._createSystem()

        if any(update.values()):
            self.energies, self.overlaps, self.vectors = None, None, None
            return True
        return False


def atom_from_config(_config):
    """Creating an Atom object from a config dictionary."""
    config = Config(_config)
    if config.nAtoms() == 1:
        return AtomOne(config)
    elif config.nAtoms() == 2:
        return AtomTwo(config)
