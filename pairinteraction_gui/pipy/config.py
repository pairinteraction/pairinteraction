"""Defining class Config.
This class is used to store all the configurations that are needed for the calculations.
"""
from hashlib import blake2b
from typing import Any
from typing import Union

import numpy as np
from pipy.misc import CustomDict

InvalidKey = object()


class Config:
    """Class for saving all the configurations that are needed for one pairinteraction calculation.
    This class is used for AtomOne as well as AtomTwo objects.

    The idea is to provide one user dictionary (dic), which is sufficient to determine all the parameters,
    but might need some calculation/tweeking to get parameters used to initialize a pairinteraction calculation.
    Therefore all parameters, which are needed for the calculation are properties of this class and
    should be accessed via the properties and not via the user dictionary!
    You can also get a complete dict of all the parameters for the calculation via the toOutDict.
    """

    def __init__(self, dic):
        """Initialize the Config class with the user dictionary dic.

        Args:
            dic: Dictionary with all the user parameters.
        """
        if not isinstance(dic, CustomDict):
            dic = CustomDict(dic, flattened=True)
        self._dic = dic  # private attribute, you should not access this directly from outside!

    def _get(self, key, default=InvalidKey):
        return self._dic.get(key, default)

    def nAtoms(self) -> int:
        """Number of involved atoms.
        Either 1 or 2, defines wether to use AtomOne or AtomTwo.
        """
        if getattr(self, "force_nAtoms", None) is not None:
            return self.force_nAtoms
        nAtoms = 1
        if self.useSameAtom() or any(k in self._dic for k in ["atom2", "species2", "atom2.species"]):
            nAtoms = 2
        return self._get("nAtoms", nAtoms)

    def useSameAtom(self) -> bool:
        """Defines if atom1 and atom2 should use the same basis.
        If True we only compute and store one single AtomOne to init AtomTwo.
        """
        return self._get("useSameAtom", self._get("samebasis", False))

    def isReal(self) -> bool:
        """Defines if we use pireal or picomplex."""
        if self.nAtoms() == 1:
            return self.isRealSingle()
        else:  # self.nAtoms() == 2
            return self.isRealSingle(1) and self.isRealSingle(2)

    def pathCache(self) -> str:
        """Path to the cache directory."""
        return self._get("pathCache", "./cache/")

    def method(self) -> str:
        """Method to be used for the calculation of the radial wavefunctions.
        Either "NUMEROV" or "WHITTAKER".
        """
        return self._get("method", "NUMEROV").upper()

    def diagonalizeThreshold(self) -> float:
        """Threshold for the diagonalization of the matrix."""
        return self._get("diagonalizeThreshold", 1e-10)

    def diamagnetism(self) -> bool:
        """Defines if we use diamagnetic corrections."""
        return self._get("diamagnetism", False)

    def toOutDict(self) -> "dict[str, Any]":
        """Get a dictionary with all the output parameters.

        Returns:
            dict: Flat dictionary with all the output parameters.
        """
        output = {
            k: getattr(self, k)()
            for k in ["nAtoms", "isReal", "method", "diagonalizeThreshold", "diamagnetism", "useSameAtom"]
        }
        output["pathCache"] = self.pathCache()
        if self.nAtoms() == 1 or self.useSameAtom():
            output.update({f"atom{getattr(self, 'iAtom', 1)}.{k}": v for k, v in self.toOutDictSingle().items()})
        else:
            for part in [1, 2]:
                output.update({f"atom{part}.{k}": v for k, v in self.toOutDictSingle(part).items()})
        if self.nAtoms() == 2:
            output.update({f"pair.{k}": v for k, v in self.toOutDictPair().items()})
        return output

    def toHash(self) -> str:
        """Return a hash of the toOutDict."""
        pairs = sorted((k, v) for k, v in self.toOutDict().items())
        h = blake2b(digest_size=10)
        h.update(str(pairs).encode("utf-8"))
        return h.hexdigest()

    def shallowCopy(self) -> "Config":
        """Just copy the Config object, but keeping the underlying user dictionary the same.
        This is useful to keep record of the original user dictionary.
        """
        return Config(self._dic)

    def deepCopy(self) -> "Config":
        """Copy the Config object and also copy the underlying user dictionary."""
        return Config(self._dic.copy())

    # Single atom parameters
    def _getSingle(self, key, default=InvalidKey, part=None):
        """Get a value from the user dictionary which corresponds to a parameter for AtomOne.
        In general options for AtomOne are stored in the user dictionary
        under the key `atom1.key`, but we also check `key1`, `atom.key` and `key` (analog for 2).
        """
        part = getattr(self, "iAtom", 1) if part is None else part
        for k in [f"atom{part}.{key}", f"{key}{part}", f"atom.{key}", key]:
            if k in self._dic:
                return self._dic[k]
        return default

    def isRealSingle(self, part=None) -> bool:
        """Wether Ey and By are zero or not."""
        isReal = self.Ey(part) == 0 and self.By(part) == 0
        return self._getSingle("isReal", isReal, part)

    def species(self, part=None) -> str:
        """Species of the atom (e.g. "Rb")."""
        return self._getSingle("species", InvalidKey, part)

    def Efield(self, part=None) -> "list[float]":
        """Electric field vector."""
        return self._getSingle("Efield", [self._getSingle("E" + x, 0, part) for x in "xyz"], part)

    def Ex(self, part=None) -> float:
        """Electric field in x direction."""
        return self.Efield(part)[0]

    def Ey(self, part=None) -> float:
        """Electric field in y direction."""
        return self.Efield(part)[1]

    def Ez(self, part=None) -> float:
        """Electric field in z direction."""
        return self.Efield(part)[2]

    def Bfield(self, part=None) -> "list[float]":
        """Magnetic field vector."""
        return self._getSingle("Bfield", [self._getSingle("B" + x, 0, part) for x in "xyz"], part)

    def Bx(self, part=None) -> float:
        """Magnetic field in x direction."""
        return self.Bfield(part)[0]

    def By(self, part=None) -> float:
        """Magnetic field in y direction."""
        return self.Bfield(part)[1]

    def Bz(self, part=None) -> float:
        """Magnetic field in z direction."""
        return self.Bfield(part)[2]

    def qnumbersSingle(self, part=None) -> "list[tuple[int, int, float, float]]":
        """Get the quantum numbers for the atom.
        Return a list of tuples (n, l, j, m).
        """
        part = getattr(self, "iAtom", 1) if part is None else part
        return self.qnumbers(part)

    def _getQnumbersSingle(self, part) -> "list[tuple[int, int, float, float]]":
        """Initialize the quantum numbers for the atom.
        Return a list of tuples (n, l, j, m).
        """
        qnumbers = self._getSingle("qnumbers", InvalidKey, part)
        if qnumbers is not InvalidKey:
            shape = np.shape(qnumbers)
            if shape == (4,):
                qnumbers = [qnumbers]
            elif len(shape) == 2 and shape[1] == 4:
                pass
            else:
                raise ValueError("Atom qnumbers must be a list of 4-tuples or a Nx4 array.")
        else:
            qnumbers = [tuple(self._getSingle(q, InvalidKey, part) for q in "nljm")]
            if any(x is InvalidKey for x in qnumbers[0]):
                raise ValueError(f"Atom qnumbers not specified for part {part}.")

        qnumbers = [(int(n), int(l), float(j), float(m)) for (n, l, j, m) in qnumbers]
        if part == 1 and self.useSameAtom():
            qnumbers += self._getQnumbersSingle(2)
        return qnumbers

    def conserveMomentaSingle(self, part=None) -> bool:
        """Wether to conserve the AtomOne momenta or not."""
        part = getattr(self, "iAtom", 1) if part is None else part
        return self.conserveMomenta(part)

    def momentaSingle(self, part=None) -> "list[float]":
        """Get the momenta, which should be conserved for the atom."""
        part = getattr(self, "iAtom", 1) if part is None else part
        return self.momenta(part)

    def restrictQnSingle(self, Q, part=None, **kwargs) -> "tuple[float, float]":
        part = getattr(self, "iAtom", 1) if part is None else part
        return self.restrictQn(Q, part, **kwargs)

    def restrictEnergySingle(self, part=None, **kwargs) -> "tuple[float, float]":
        return self.restrictQnSingle("Energy", part, **kwargs)

    def getEnergiesSingle(self, part=None, pi=None) -> "list[float]":
        """Get the energies of the atom states in qnumbersSingle without any E/B fields."""
        part = getattr(self, "iAtom", 1) if part is None else part
        if getattr(self, f"_energiesSingle{part}", None) is None:
            setattr(self, f"_energiesSingle{part}", self._getSingle("_energies", None, part))
        if getattr(self, f"_energiesSingle{part}", None) is None:
            if pi is None:
                import pipy

                pi = pipy.pireal
            states = [pi.StateOne(self.species(), *qn) for qn in self.qnumbersSingle(part)]
            energies = [s.getEnergy() for s in states]
            setattr(self, f"_energiesSingle{part}", energies)
        return getattr(self, f"_energiesSingle{part}")

    def toOutDictSingle(self, part=None) -> "dict[str, Any]":
        """Get a dictionary with all the output parameters for one AtomOne."""
        output = {
            "species": self.species(part),
            "momenta": self.momentaSingle(part),
        }
        output["_energies"] = getattr(self, f"_energiesSingle{part}", None)
        for Q in ["Energy", "N", "L", "J", "M"]:
            minQ, maxQ = self.restrictQnSingle(Q, part)
            output.update({f"min{Q}": minQ, f"max{Q}": maxQ})
        for x, v in zip("xyz", self.Efield(part)):
            output[f"E{x}"] = v
        for x, v in zip("xyz", self.Bfield(part)):
            output[f"B{x}"] = v
        return output

    # Pair atom parameters
    def _getPair(self, key, default=InvalidKey):
        """Get a value from the user dictionary which corresponds to a parameter for AtomTwo.
        In general options for AtomTwo are stored in the user dictionary
        under the key `pair.key` but we also check `key`.
        """
        for k in [f"pair.{key}", key]:
            if k in self._dic:
                return self._dic[k]
        return default

    def order(self) -> int:
        """Get the order, up to which the multipole expansion should be calculated."""
        return self._getPair("order", self._getPair("exponent", 3))

    def minimalNorm(self) -> float:
        """Get the minimal norm for states to be included in the basis of SystemTwo."""
        return self._getPair("minimalNorm", 5e-2)

    def angle(self) -> float:
        """Get the angle between the two atoms in radian."""
        return self._getPair("angle", self._getPair("theta", 0))

    def distance(self) -> float:
        r"""Get the distance between the two atoms \mu m."""
        return self._getPair("distance", self._getPair("R", np.inf))

    def symmetry(self, sym) -> Union[None, str]:
        """Get a symmetry of the system.
        Returns "ODD", "EVEN" or None.
        """
        assert sym in ["inversion", "permutation", "reflection"]
        val = self._getPair(sym, self._getPair(sym[:3], None))
        val = None if val is None else val.upper()
        dic = {
            None: [None, "NONE"],
            "ODD": ["ODD", "O"],
            "EVEN": ["EVEN", "E"],
        }
        for k, v in dic.items():
            if val in v:
                return k
        raise ValueError(f"Invalid symmetry value for {sym}: {val}")

    def inversion(self) -> Union[None, str]:
        return self.symmetry("inversion")

    def permutation(self) -> Union[None, str]:
        return self.symmetry("permutation")

    def reflection(self) -> Union[None, str]:
        return self.symmetry("reflection")

    def useQnumbersSingleAsPair(self) -> bool:
        """If True you should only provide in total two qnumbersSingle for the two atoms.
        The combination of the two qnumbersSingle will then be used as the qnumbers for the pair.
        """
        return self._getPair("useQnumbersSingleAsPair", False)

    def qnumbersPair(self) -> "list[tuple[int, int, float, float, int, int, float, float]]":
        """Get the quantum numbers for the pair.
        Return a list of tuples (n1, l1, j1, m1, n2, l2, j2, m2).
        """
        return self.qnumbers("pair")

    def _getQnumbersPair(self) -> "list[tuple[int, int, float, float, int, int, float, float]]":
        """Initialize the quantum numbers for the pair.
        Return a list of tuples (n1, l1, j1, m1, n2, l2, j2, m2).
        """
        qnumbers = self._getPair("qnumbers", InvalidKey)
        if qnumbers is not InvalidKey:
            shape = np.shape(qnumbers)
            if shape == (8,):
                qnumbers = [qnumbers]
            elif len(shape) == 2 and shape[1] == 8:
                pass
            elif len(shape) == 3 and shape[1:] == (4, 2):
                # reshape [(N1, N2), (L1, L2), (J1, J2), (M1, M2)] to (N1, L1, J1, M1, N2, L2, J2, M2)
                qnumbers = np.reshape(qnumbers, (shape[0], 8), order="F")
            elif len(shape) == 3 and shape[1:] == (2, 4):
                # reshape [(N1, L1, J1, M1), (N2, L2, J2, M2)] to (N1, L1, J1, M1, N2, L2, J2, M2)
                qnumbers = np.reshape(qnumbers, (shape[0], 8), order="C")
            else:
                raise ValueError("Atom qnumbers must be a list of 4-tuples or a Nx4 array.")
        elif self.useQnumbersSingleAsPair():
            qnsSingle = self.qnumbersSingle(1)
            if not self.useSameAtom():
                qnsSingle += self.qnumbersSingle(2)
            if len(qnsSingle) != 2:
                raise ValueError("When useQnumbersSingleAsPair there must be exactly two qnumbersSingle.")
            qnumbers = [(*qnsSingle[0], *qnsSingle[1])]
        else:
            qnumbers = []
            for i, qn1 in enumerate(self.qnumbersSingle(1)):
                for qn2 in self.qnumbersSingle(2)[i * self.useSameAtom() :]:
                    qnumbers.append((*qn1, *qn2))

        qnumbers = [
            (int(n1), int(l1), float(j1), float(m1), int(n2), int(l2), float(j2), float(m2))
            for (n1, l1, j1, m1, n2, l2, j2, m2) in qnumbers
        ]
        return qnumbers

    def conserveMomentaPair(self) -> bool:
        """Wether to conserve the AtomTwo momenta or not."""
        return self.conserveMomenta("pair")

    def momentaPair(self) -> "list[float]":
        """Get the momenta, which should be conserved for the pair."""
        return self.momenta("pair")

    def restrictQnPair(self, Q, **kwargs) -> "tuple[float, float]":
        """Allowed formatting in config:
        - dQ, deltaQ, dQPair, deltaQPair
        - minQ, minQPair, maxQ, maxQPair
        """
        return self.restrictQn(Q, "pair", **kwargs)

    def restrictEnergyPair(self, **kwargs) -> "tuple[float, float]":
        return self.restrictQnPair("Energy", **kwargs)

    def getEnergiesPair(self, atom1=None, atom2=None) -> "list[float]":
        if getattr(self, "_energiesPair", None) is None:
            self._energiesPair = self._getPair("_energies", None)
        if getattr(self, "_energiesPair", None) is None:
            if atom1 is None or atom2 is None:
                raise NotImplementedError("For now atom1 and atom2 must be given to calculate Pair energies.")
                # TODO: implement this via import AtomOne ...
            qns1 = np.array(self.qnumbersPair())[:, :4]
            qns2 = np.array(self.qnumbersPair())[:, 4:]
            E1 = atom1.getEnergiesStates(qns1)
            E2 = atom2.getEnergiesStates(qns2)
            self._energiesPair = list(E1 + E2)
        return self._energiesPair

    def toOutDictPair(self) -> "dict[str, Any]":
        """Get a dictionary with all the output parameters for one AtomOne."""
        output = {
            "order": self.order(),
            "minimalNorm": self.minimalNorm(),
            "momenta": self.momentaPair(),
            "angle": self.angle(),
            "distance": self.distance(),
        }
        output["_energies"] = getattr(self, "_energiesPair", None)
        for Q in ["Energy"]:  # ["Energy", "N", "L", "J", "M"]:
            minQ, maxQ = self.restrictQnPair(Q)
            output.update({f"min{Q}": minQ, f"max{Q}": maxQ})
        for sym in ["inversion", "permutation", "reflection"]:
            output[sym] = self.symmetry(sym)
        return output

    # Single and Pair atom properties
    # Many parameters are defined for both single and pair atoms
    # To avoid code duplication we define the following functions
    # Here the parameter part defines if we are looking for a single (and which atom) or pair parameter
    # The parameter part can be one of the following:
    # - None: We are looking for a single parameter and just a single atom is involved (will default to 1)
    # - 1 or 2: We are looking for a single parameter and the atom with the given index is involved
    # - "pair": We are looking for a pair parameter
    def _getPart(self, key, default=InvalidKey, part=None):
        """Wrapper around _getSingle and _getPair."""
        if str(part).lower() == "pair":
            return self._getPair(key, default)
        return self._getSingle(key, default, part)

    def qnumbers(self, part) -> "list[tuple]":
        """Get the quantum numbers for the part.
        Return a list of tuples.
        """
        name = "Pair" if str(part).lower() == "pair" else f"Single{part}"
        if not hasattr(self, f"_qnumbers{name}"):
            setattr(self, f"_qnumbers{name}", self._getQnumbers(part))
        return getattr(self, f"_qnumbers{name}")

    def _getQnumbers(self, part) -> "list[tuple]":
        """Wrapper around _getQnumbersSingle and _getQnumbersPair."""
        if str(part).lower() == "pair":
            return self._getQnumbersPair()
        return self._getQnumbersSingle(part)

    def conserveMomenta(self, part) -> bool:
        """Wether to conserve the momenta or not."""
        conserveMomenta = self._getPart("momenta", None, part) is not None
        return self._getPart("conserveMomenta", conserveMomenta, part)

    def momenta(self, part) -> "list":
        """Get the momenta, which should be conserved for the atom."""
        ispair = str(part).lower() == "pair"
        momenta = self._getPart("momenta", None, part)
        if momenta is not None:
            if isinstance(momenta, (int, float)):
                momenta = [momenta]
        elif self.conserveMomenta(part):
            qnumbers = np.array(self.qnumbers(part))
            Ms = qnumbers[:, 3] + qnumbers[:, 7] if ispair else qnumbers[:, 3]
            if not np.all(Ms == Ms[0]):
                raise ValueError("Not all atomstates have the same M, cannot conserve momenta.")
            momenta = [Ms[0]]
        else:
            return None
        return [int(m) if ispair else float(m) for m in momenta]

    def restrictQn(self, Q, part=None, **kwargs) -> "tuple[Union[float, int], Union[float, int]]":
        """Allowed formatting in config:
        - dQ, deltaQ, dQSingle, deltaQSingle
        - minQ, minQSingle, maxQ, maxQSingle
        """
        Q = Q.capitalize() if Q.capitalize() != "E" else "Energy"
        ispair = str(part).lower() == "pair"
        single_pair = "Pair" if ispair else "Single"
        key_list = [f"d{Q}", f"delta{Q}", f"d{Q}{single_pair}", f"delta{Q}{single_pair}"]
        if Q == "Energy":
            key_list += [f"deltaE{single_pair}"]
        for key in key_list:
            dQ = self._getPart(key, InvalidKey, part)
            if dQ is not InvalidKey:
                break
        dQ = None if (dQ is InvalidKey or dQ is None or dQ < 0) else dQ
        minQ = self._getPart(f"min{Q}", None, part)
        maxQ = self._getPart(f"max{Q}", None, part)

        if minQ is not None and maxQ is not None:
            pass
        elif dQ is not None:
            if Q == "Energy":
                values = self.getEnergies(part, **kwargs)
            else:
                index = ["N", "L", "J", "M"].index(Q)
                qnumbers = np.array(self.qnumbers(part))
                values = qnumbers[:, index]
                values = qnumbers[:, index] + qnumbers[:, index + 4] if ispair else qnumbers[:, index]
            minQ, maxQ = min(values) - dQ, max(values) + dQ
        else:
            return None, None
        _type = int if Q in "NL" or (ispair and Q in "JM") else float
        return _type(minQ), _type(maxQ)

    def getEnergies(self, part, **kwargs) -> "list[float]":
        """Wrapper around getEnergiesSingle and getEnergiesPair."""
        if str(part).lower() == "pair":
            return self.getEnergiesPair(**kwargs)
        return self.getEnergiesSingle(part, **kwargs)

    # set stuff
    def setEfield(self, value, part):
        self._dic[f"atom{part}.Efield"] = value

    def setBfield(self, value, part):
        self._dic[f"atom{part}.Bfield"] = value

    def setSymmetry(self, sym, val):
        self._dic["pair." + sym] = val

    def setAngle(self, val):
        self._dic["pair.angle"] = val

    def setDistance(self, val):
        self._dic["pair.distance"] = val
