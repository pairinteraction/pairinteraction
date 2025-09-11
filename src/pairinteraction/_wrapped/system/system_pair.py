# SPDX-FileCopyrightText: 2024 PairInteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

import logging
from typing import TYPE_CHECKING, Any, ClassVar, Optional, TypeVar, Union, overload

import numpy as np

from pairinteraction import _backend
from pairinteraction._wrapped.basis.basis_pair import BasisPairComplex, BasisPairReal
from pairinteraction._wrapped.system.green_tensor import GreenTensor
from pairinteraction._wrapped.system.system import SystemBase
from pairinteraction.units import QuantityScalar

if TYPE_CHECKING:
    from typing_extensions import Self

    from pairinteraction._wrapped.basis.basis_pair import BasisPair
    from pairinteraction._wrapped.ket.ket_atom import (
        KetAtom,  # noqa: F401  # needed for sphinx to recognize KetAtomTuple
    )
    from pairinteraction._wrapped.ket.ket_pair import KetAtomTuple
    from pairinteraction.units import (
        ArrayLike,
        PintArray,  # noqa: F401  # needed for sphinx to recognize PintArrayLike
        PintArrayLike,
        PintFloat,
    )

BasisType = TypeVar("BasisType", bound="BasisPair[Any, Any]", covariant=True)
UnionCPPSystemPair = Union[_backend.SystemPairReal, _backend.SystemPairComplex]
UnionTypeCPPSystemPair = Union[type[_backend.SystemPairReal], type[_backend.SystemPairComplex]]

logger = logging.getLogger(__name__)


class SystemPair(SystemBase[BasisType]):
    """System of a pair of atoms.

    Use the given BasisPair object to create the system object.
    You can set the distance (vector) between the atoms afterwards via the corresponding methods.

    Examples:
        >>> import pairinteraction.real as pi
        >>> ket = pi.KetAtom("Rb", n=60, l=0, m=0.5)
        >>> basis = pi.BasisAtom("Rb", n=(58, 63), l=(0, 3))
        >>> system = pi.SystemAtom(basis).set_electric_field([0.1, 0, 0.1], unit="V/cm").diagonalize()
        >>> system = pi.SystemAtom(basis).set_magnetic_field([0, 0, 1], unit="G").diagonalize()
        >>> pair_energy = 2 * system.get_corresponding_energy(ket, unit="GHz")
        >>> pair_basis = pi.BasisPair(
        ...     [system, system],
        ...     energy=(pair_energy - 3, pair_energy + 3),
        ...     energy_unit="GHz",
        ... )
        >>> pair_system = pi.SystemPair(pair_basis).set_distance(5, unit="micrometer").set_interaction_order(3)
        >>> print(pair_system)
        SystemPair(BasisPair(|Rb:59,S_1/2,-1/2; Rb:61,S_1/2,-1/2⟩ ... |~Rb:58,F_7/2,7/2; Rb:59,S_1/2,1/2⟩), is_diagonal=False)
        >>> pair_system = pair_system.diagonalize()
        >>> eigenenergies = pair_system.get_eigenenergies(unit="GHz")
        >>> print(f"{eigenenergies[0] - pair_energy:.5f}")
        -2.18394

    """  # noqa: E501

    _cpp: UnionCPPSystemPair
    _cpp_type: ClassVar[UnionTypeCPPSystemPair]

    def __init__(self, basis: BasisType) -> None:
        """Create a system object for a pair of atoms.

        Args:
            basis: The :class:`pairinteraction.real.BasisPair` object that describes the basis of the system.

        """
        super().__init__(basis)
        self._distance_vector_au = [0, 0, np.inf]
        self._interaction_order = 3

    def get_eigenbasis(self) -> BasisType:
        """Get the eigenbasis of the system.

        This method retrieves the eigenbasis of the system, which is the basis in which the Hamiltonian is diagonal.

        Returns:
            A BasisPair object representing the eigenbasis of the system.

        """
        eigenbasis = super().get_eigenbasis()
        eigenbasis.system_atoms = self.basis.system_atoms
        return eigenbasis

    def _update_basis(self) -> None:
        system_atoms = self.basis.system_atoms
        super()._update_basis()
        self._basis.system_atoms = system_atoms

    def set_interaction_order(self: "Self", order: int) -> "Self":
        """Set the interaction order of the pair system.

        Args:
            order: The interaction order to set for the pair system.
                The order must be 3, 4, or 5.

        """
        self._interaction_order = order
        self._cpp.set_interaction_order(order)
        return self

    def set_distance(
        self: "Self",
        distance: Union[float, "PintFloat"],
        angle_degree: float = 0,
        unit: Optional[str] = None,
    ) -> "Self":
        """Set the distance between the atoms using the specified distance and angle.

        Args:
            distance: The distance to set between the atoms in the given unit.
            angle_degree: The angle between the distance vector and the z-axis in degrees.
                90 degrees corresponds to the x-axis.
                Defaults to 0, which corresponds to the z-axis.
            unit: The unit of the distance, e.g. "micrometer".
                Default None expects a `pint.Quantity`.

        """
        distance_vector = [np.sin(np.deg2rad(angle_degree)) * distance, 0, np.cos(np.deg2rad(angle_degree)) * distance]
        return self.set_distance_vector(distance_vector, unit)

    def set_distance_vector(
        self: "Self",
        distance: Union["ArrayLike", "PintArrayLike"],
        unit: Optional[str] = None,
    ) -> "Self":
        """Set the distance vector between the atoms.

        Args:
            distance: The distance vector to set between the atoms in the given unit.
            unit: The unit of the distance, e.g. "micrometer".
                Default None expects a `pint.Quantity`.

        """
        distance_au = [QuantityScalar.convert_user_to_au(v, unit, "distance") for v in distance]
        self._cpp.set_distance_vector(distance_au)
        self._distance_vector_au = distance_au
        return self

    @overload
    def get_distance_vector(self, unit: None = None) -> list["PintFloat"]: ...

    @overload
    def get_distance_vector(self, unit: str) -> list[float]: ...

    def get_distance_vector(self, unit: Optional[str] = None) -> Union[list[float], list["PintFloat"]]:
        return [QuantityScalar.convert_au_to_user(d, "distance", unit) for d in self._distance_vector_au]  # type: ignore [return-value]

    @overload
    def get_distance(self, unit: None = None) -> "PintFloat": ...

    @overload
    def get_distance(self, unit: str) -> float: ...

    def get_distance(self, unit: Optional[str] = None) -> Union[float, "PintFloat"]:
        distance = np.linalg.norm(self._distance_vector_au)
        return QuantityScalar.convert_au_to_user(float(distance), "distance", unit)

    def set_green_tensor(self, green_tensor: GreenTensor) -> "Self":
        """Set the Green tensor for the pair system.

        Args:
            green_tensor: The Green tensor to set for the system.

        """
        if isinstance(self._cpp, _backend.SystemPairReal) and isinstance(green_tensor._cpp, _backend.GreenTensorReal):
            self._cpp.set_green_tensor(green_tensor._cpp)
        elif isinstance(self._cpp, _backend.SystemPairComplex) and isinstance(
            green_tensor._cpp, _backend.GreenTensorComplex
        ):
            self._cpp.set_green_tensor(green_tensor._cpp)
        else:
            raise TypeError(f"Incompatible types: {type(self)=}; {type(green_tensor)=}")

        return self

    @overload
    def get_corresponding_energy(self: "Self", ket_tuple: "KetAtomTuple", unit: None = None) -> "PintFloat": ...

    @overload
    def get_corresponding_energy(self: "Self", ket_tuple: "KetAtomTuple", unit: str) -> float: ...

    def get_corresponding_energy(
        self: "Self", ket_tuple: "KetAtomTuple", unit: Optional[str] = None
    ) -> Union[float, "PintFloat"]:
        overlaps = self.get_eigenbasis().get_overlaps(ket_tuple)
        idx = np.argmax(overlaps)
        if overlaps[idx] <= 0.5:
            logger.warning(
                "The provided ket states does not correspond to an eigenstate of the system in a unique way."
            )
        return self.get_eigenenergies(unit=unit)[idx]  # type: ignore [index,no-any-return] # PintArray does not know it can be indexed


class SystemPairReal(SystemPair[BasisPairReal]):
    _cpp: _backend.SystemPairReal
    _cpp_type = _backend.SystemPairReal
    _TypeBasis = BasisPairReal


class SystemPairComplex(SystemPair[BasisPairComplex]):
    _cpp: _backend.SystemPairComplex
    _cpp_type = _backend.SystemPairComplex
    _TypeBasis = BasisPairComplex
