# SPDX-FileCopyrightText: 2025 PairInteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later
from __future__ import annotations

from typing import TYPE_CHECKING, overload

from pairinteraction.perturbative.effective_system_pair import EffectiveSystemPair, EffectiveSystemPairReal
from pairinteraction.units import QuantityScalar

if TYPE_CHECKING:
    from pairinteraction.ket import KetAtom
    from pairinteraction.units import PintFloat


class C6(EffectiveSystemPair):
    """Class for calculating the C6 coefficient for a two-atom state ``|ket1, ket2>``.

    Given two `KetAtom` objects ket1 and ket2,
    this class computes the C6 coefficient for the pair state ``|ket1, ket2>``.
    This class also allows to set magnetic and electric fields similar to the `SystemAtom` class,
    as well as the angle between the two atoms like in the `SystemPair` class.

    Note, that ket1 and ket2 must be either the same state or states of different species.
    If you want to calculate the C6 coefficient for two different states of the same species,
    we recommend using the `EffectiveSystemPair` class with the ket_tuples subspace [(ket1, ket2), (ket2, ket1)].



    Examples:
        >>> import pairinteraction as pi
        >>> ket = pi.KetAtom("Rb", n=60, l=0, j=0.5, m=0.5)
        >>> c6_obj = pi.C6(ket, ket)
        >>> c6 = c6_obj.get(unit="planck_constant * GHz * micrometer^6")
        >>> print(f"{c6:.1f}")
        138.9

    """

    def __init__(self, ket1: KetAtom, ket2: KetAtom) -> None:
        if ket1.species == ket2.species and ket1 != ket2:
            raise ValueError(
                "If you want to calculate 2nd order perturbations of two different states a and b "
                "(of the same species), "
                "please use the EffectiveSystemPair([(a,b), (b,a)]) class."
            )

        super().__init__([(ket1, ket2)])
        self.set_perturbation_order(2)

        # Set some default distance. The exact value does not matter for the C6 value,
        # but too small distances result in warnings about the basisvec admixtures
        self.set_distance(100, unit="micrometer")

    @overload
    def get(self, unit: None = None) -> PintFloat: ...

    @overload
    def get(self, unit: str) -> float: ...

    def get(self, unit: str | None = None) -> float | PintFloat:
        """Get the C6 coefficient of the pair interaction between the specified ket1 and ket2.

        Args:
            unit: The unit in which to return the C6 coefficient.
                If None, returns a pint object.

        """
        h_eff_pint = self.get_effective_hamiltonian(return_order=2)
        distance = self.system_pair.get_distance()
        c6_pint = h_eff_pint[0, 0] * distance**6  # type: ignore [index]  # pint does not know it can be indexed
        return QuantityScalar.convert_pint_to_user(c6_pint, "c6", unit)


class C6Real(C6, EffectiveSystemPairReal):
    pass
