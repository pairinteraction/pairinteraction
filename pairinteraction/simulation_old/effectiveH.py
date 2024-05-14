"""Defining class EffectiveH
"""

from typing import Optional, Union

import numpy as np

from pairinteraction import pireal
from pairinteraction.model.model import ModelSimulation
from pairinteraction.simulation_old.atom import AtomTwo
from pairinteraction.simulation_old.exceptions import CppDeleted, CppObjectAlreadyDeleted


class EffectiveH:
    def __init__(self, model: Union[ModelSimulation, dict], options: Optional[dict] = None):
        if not isinstance(model, ModelSimulation):
            model = ModelSimulation.model_validate(model)
        self.model = model
        self.options = {} if options is None else options

        # this atomtwo is only used as shorthand access to it's atom1 and atom2,
        # since we did not choose any subspace yet, we never should build the atomtwo system!
        self.atom = AtomTwo(self.model)
        self.atom._system = "Do not access this!"

        self.hamiltonians = {}
        self.basis = None

    def __repr__(self):
        return f"{self.__class__.__name__}(...)"

    def run(self):
        """Run all relevant calculations and return the constructed Hamiltonians.

        Returns
        -------
        Hamiltonians: dict
            A dictionary containing the Hamiltonians, with keys "unperturbed", "perturbed", "effective", "onsite".
        """
        self.initSubspaces()
        self.calcSubspaces()
        self.constructTotalHamiltonians()
        self.constructOnsite()
        return self.hamiltonians

    def initSubspaces(self):
        """Initializes self.subspacesIds.

        subspacesIds: list of list of tuples (id1, id2),
        where id1 (id2) are the index of the state to use from the atom1 (atom2) basis
        """
        if self.options.get("subspacesIds", None) is not None:
            self.subspacesIds = self.options["subspacesIds"]
            self.Nsubspaces = len(self.subspacesIds)
            self.dim = len([idSingle for ids in self.subspacesIds for idSingle in ids])
            return

        # else: init all possible subspaces from atom1 and atom2 qnumbers,
        # which twoatom states are in the twoatom dEnergy range
        atom, model = self.atom, self.model
        dEnergy = model.interactions.delta_energy
        if dEnergy is None:
            raise ValueError("dEnergy not set, cannot calculate subspaces!")

        Es = {
            1: [atom.atom1.getStateEnergy(s) for s in model.atom1.states_of_interest],
            2: [atom.atom2.getStateEnergy(s) for s in model.atom2.states_of_interest],
        }

        EPair = [E1 + E2 for E1 in Es[1] for E2 in Es[2]]
        idSingle = [(id1, id2) for id1 in range(len(Es[1])) for id2 in range(len(Es[2]))]

        subspaces = []
        oldEnergy = np.inf
        for idPair in np.argsort(EPair):
            if np.abs(oldEnergy - EPair[idPair]) > dEnergy:
                subspaces.append([])
            subspaces[-1].append(idSingle[idPair])
            oldEnergy = EPair[idPair]

        self.Nsubspaces = len(subspaces)
        self.subspacesIds = subspaces
        self.dim = len(idSingle)

    def calcSubspaces(self):
        self.effectiveSubspaces = []
        for ids in self.subspacesIds:
            model = self.model.model_dump()
            model["interactions"]["combined_states_of_interest"] = [{"atom1": i1, "atom2": i2} for (i1, i2) in ids]
            effS = EffectiveSubspace(model)
            effS.run()
            self.effectiveSubspaces.append(effS)

    def constructTotalHamiltonians(self):
        for pert in ["unperturbed", "perturbed"]:
            dtype = np.common_type(*[effS.hamiltonians[pert] for effS in self.effectiveSubspaces])
            H = np.zeros((self.dim, self.dim), dtype=dtype)

            # Create total Hamiltonian
            end = 0
            for effS in self.effectiveSubspaces:
                H_subspace = effS.hamiltonians[pert]
                start = end
                end = start + H_subspace.shape[0]
                H[start:end, start:end] = H_subspace

            # Reorder the basis of the total Hamiltonian
            oldBasis = [idSingle for ids in self.subspacesIds for idSingle in ids]
            dim2 = len(self.model.atom2.states_of_interest)
            basisMapping = np.argsort([id1 * dim2 + id2 for (id1, id2) in oldBasis])
            self.hamiltonians[pert] = H[basisMapping[:, np.newaxis], basisMapping]

        self.hamiltonians["effective"] = self.hamiltonians["perturbed"] - self.hamiltonians["unperturbed"]
        self.basisIds = [oldBasis[i] for i in basisMapping]
        self.basis = [
            (self.model.atom1.states_of_interest[id1], self.model.atom2.states_of_interest[id2])
            for (id1, id2) in self.basisIds
        ]

    def constructOnsite(self):
        """Construct onsite energies"""
        assert self.model.atom1 == self.model.atom2
        Es = [self.atom.atom1.getStateEnergy(s) for s in self.model.atom1.states_of_interest]
        self.hamiltonians["onsite"] = np.diag(Es)

    def deleteCppObjects(self):
        """Delete all cpp objects, usually you want to call this, before saving the results
        to a pkl file to make the size much smaller and prevent weird bugs due to the cpp objects.
        """
        # Never delete _cache without also deleting the _system!
        # Also first delete _system, then _cache
        self.atom.deleteCppObjects()
        for effS in self.effectiveSubspaces:
            effS.deleteCppObjects()


class EffectiveSubspace:
    """Class describing one subspace of the EffectiveH class."""

    def __init__(self, model: Union[ModelSimulation, dict]):
        if not isinstance(model, ModelSimulation):
            model = ModelSimulation.model_validate(model)
        self.model = model
        self.configInt = {
            "distance": model.interactions.distance.get_value(),
            "angle": model.interactions.angle.get_value(),
        }
        model.interactions.distance = None
        model.interactions.angle = None

        self.atom = AtomTwo(model)
        self.hamiltonians = {}

    def run(self):
        self._buildSystemUnperturbed()
        self._buildSystemPerturbed()

    @property
    def systemUnperturbed(self):
        if getattr(self, "_systemUnperturbed", None) is None:
            self._buildSystemUnperturbed()
        elif self._systemUnperturbed is CppDeleted:
            raise CppObjectAlreadyDeleted(self.__class__.__name__ + ".systemUnperturbed")
        return self._systemUnperturbed

    @property
    def phaseCorrection(self):
        if getattr(self, "_phaseCorrection", None) is None:
            self._buildSystemUnperturbed()
        return self._phaseCorrection

    @property
    def systemPerturbed(self):
        if getattr(self, "_systemPerturbed", None) is None:
            self._buildSystemPerturbed()
        elif self._systemPerturbed is CppDeleted:
            raise CppObjectAlreadyDeleted(self.__class__.__name__ + ".systemPerturbed")
        return self._systemPerturbed

    def _buildSystemUnperturbed(self):
        """Construct unperturbed system and Hamiltonian"""
        # get the unperturbed (= non interacting) system, the two body Hamiltonian
        # (with shape (# basis states, # basis states)) is diagonal
        # (Note the single atom eigenstates are used as basis)
        system = self.atom.copySystem()  # system.getBasisvectors().shape = (# involved states, # basis states)
        # now we constrain the basis states to only consist of the wanted states
        subspace = [
            pireal.StateTwo(cs["atom1"].state, cs["atom2"].state)
            for cs in self.model.interactions.combined_states_of_interest
        ]
        system.constrainBasisvectors(
            system.getBasisvectorIndex(subspace)
        )  # system.getBasisvectors().shape = (# involved states, # wanted states)

        # Unfortunately, when building the system, the basis states can have
        # a plus or minus sign in front of the coefficient with the greatest overlap
        # Note, the basis states should represent the wanted states (self.subspacesIds[indS]),
        # thus to get the correct effective Hamiltonian for those states,
        # we want a plus sign in front of that coefficient,
        # hence we should correct this sign, which is done in _getPhaseCorrection
        # Also note, that we have to do this here, before unitarizing the system,
        # since unitarizing the system will brute force set the getBasisvectors to diag([1, 1, ..., 1]),
        # independent of the discussed sign.
        # Note 3, for the unpertubed case, this sign/phase correction is not needed,
        # since the unperturbed Hamiltonian is diagonal anyway, however it is needed for the perturbed case.
        self._phaseCorrection = self._getPhaseCorrection(system)

        # Finally we unitarize the system (i.e. create as many new involved states as there are basis states
        # (=wanted states in this case) by defining them via their overlap with the old involved states
        # and label them by some generated hash, see pairinteraction SystemBase.cpp:unitarize)
        # this is needed for SchriefferWolffTransformation

        system.unitarize()
        # system.getBasisvectors().shape = (# wanted states, # wanted states) with simply 1 on the diagonal

        self._systemUnperturbed = system

        H = system.getHamiltonian().todense()
        H = np.real_if_close(self.phaseCorrection.conjugate() @ H @ self.phaseCorrection, tol=1e-12)
        self.hamiltonians["unperturbed"] = H

    def _buildSystemPerturbed(self):
        """Construct perturbed system and Hamiltonian"""
        # build the perturbed (= interacting) system, the two body Hamiltonian is not diagonal anymore,
        # but we still use the single atom eigenstates as basis
        system = self.atom.copySystem()
        if self.configInt["distance"] is None:
            raise ValueError("Distance not set, cannot calculate perturbed system!")
        system.setDistance(self.configInt["distance"])
        if self.configInt["angle"] is not None:
            system.setAngle(self.configInt["angle"])
        # system.getBasisvectors().shape = (# involved states, # basis states)

        # Since we want to take higher order interactions into account, we cannot simply constrain the basis,
        # but we are applying the SchriefferWolffTransformation
        # (actually we should call it Direct Rotation) to the system.
        # However, for this we first need to unitarize the system.
        # This is similar to the unperturbed case, but now we need all basis states, and not just the wanted states.
        system.unitarize()

        # Finally, apply the SchriefferWolffTransformation, to get the effective Hamiltonian
        # applySchriefferWolffTransformation will also diagonalize the system at the start (with threshold 0)
        system.applySchriefferWolffTransformation(self.systemUnperturbed)

        # phase correction same as from unperturbed case, since we are "projecting"
        # to the same basis via the SchriefferWolffTransformation/Direct Rotation

        self._systemPerturbed = system

        H = system.getHamiltonian().todense()
        H = np.real_if_close(self.phaseCorrection.conjugate() @ H @ self.phaseCorrection, tol=1e-12)
        self.hamiltonians["perturbed"] = H

    def _getPhaseCorrection(self, system):
        r"""Building the basis for the two atom Hamiltonian with electric and magnetic fields
        can introduce sign changes (or complex phases if Ey, By != 0) in the coefficients of the new basis
        with the greatest overlap of the old basis.
        To correct for this signs this function looks if/which signs (phase)
        did change and returns a diagonal matrix with +/-1 (exp(i*phase)) to correct this sign (phase) error.

        In more detail this means:
        If you start with a state |ac>=|a,c> (atom1 in state |a>, atom2 in state |b>)
        and apply E and B fields you will in general end up with a state:
        :math:`|\tilde{ac}> = \sum_i c_i |state_i>`
        where for "small fields" the coefficent :math:`\abs(c_{ab}) \approx 1`.
        But there can be an internal phase change such that :math:`c_{ab} \approx exp(i*phase)`
        this function checks for these phase changes and returns a diagonal matrix with exp(i*phase)
        entries dependent on this phase, to correct for this phase change.
        This is important, because it can effect the sign of the offdiagonal elements in the effective Hamiltonian.
        """
        # TODO save also the overlap value somewhere, that could be interesting for looking,
        # wheter perturbatively calculating the effective Hamiltonian is a good idea
        basisvectors = system.getBasisvectors()
        phaseCorrection = np.identity(basisvectors.shape[1], dtype=basisvectors.dtype)
        for i in range(basisvectors.shape[1]):
            ind = np.argmax(np.abs(basisvectors[:, i]))
            val = basisvectors[ind, i]
            if np.abs(val) < 0.55:
                raise ValueError(
                    "_getPhaseCorrection: Initial state not unambiguous anymore."
                    "This should not happen when calculating the effective Hamiltonian perturbatively."
                )
            if np.imag(val) == 0:
                phaseCorrection[i, i] = 1 if val > 0 else -1
            else:
                phaseCorrection[i, i] = np.exp(-np.angle(val) * 1j)
        return phaseCorrection

    def deleteCppObjects(self):
        """Delete all cpp objects, usually you want to call this, before saving the results
        to a pkl file to make the size much smaller and prevent weird bugs due to the cpp objects.
        """
        self.atom.deleteCppObjects()

        for k in ["_systemUnperturbed", "_systemPerturbed"]:
            delattr(self, k)
            setattr(self, k, CppDeleted)
