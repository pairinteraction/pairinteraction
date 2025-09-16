# SPDX-FileCopyrightText: 2025 PairInteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

from typing import TYPE_CHECKING, Optional, Union, overload

from pairinteraction.perturbative.effective_system_pair import EffectiveSystemPair, EffectiveSystemPairReal
from pairinteraction.units import QuantityScalar

if TYPE_CHECKING:
    from pairinteraction.ket import KetAtom
    from pairinteraction.units import PintFloat


class C3(EffectiveSystemPair):
    """Class for calculating the C3 coefficient between two states.

    Given two `KetAtom` objects ket1 and ket2,
    this class computes the C3 coefficient between ``|ket1, ket2>`` and ``|ket2, ket1>``.
    This class also allows to set magnetic and electric fields similar to the `SystemAtom` class,
    as well as the angle between the two atoms like in the `SystemPair` class.

    Examples:
        >>> import pairinteraction as pi
        >>> ket1 = pi.KetAtom("Rb", n=60, l=0, j=0.5, m=0.5)
        >>> ket2 = pi.KetAtom("Rb", n=61, l=1, j=1.5, m=0.5)
        >>> c3_obj = pi.C3(ket1, ket2)
        >>> c3 = c3_obj.get(unit="planck_constant * MHz * micrometer^3")
        >>> print(f"{c3:.2f}")
        -93.29

    """

    def __init__(self, ket1: "KetAtom", ket2: "KetAtom") -> None:
        super().__init__([(ket1, ket2), (ket2, ket1)])
        self.set_perturbation_order(1)

        # Set some default distance. The exact value does not matter for the C3 value
        self.set_distance(100, unit="micrometer")

    @overload
    def get(self, unit: None = None) -> "PintFloat": ...

    @overload
    def get(self, unit: str) -> float: ...

    def get(self, unit: Optional[str] = None) -> Union[float, "PintFloat"]:
        """Get the C3 coefficient of the pair interaction between the specified ket1 and ket2.

        Args:
            unit: The unit in which to return the C3 coefficient.
                If None, returns a pint object.

        """
        h_eff_pint = self.get_effective_hamiltonian(return_order=1)
        distance = self.system_pair.get_distance()
        c3_pint = h_eff_pint[1, 0] * distance**3  # type: ignore [index]  # pint does not know it can be indexed
        return QuantityScalar.convert_pint_to_user(c3_pint, "c3", unit)


class C3Real(C3, EffectiveSystemPairReal):
    pass
