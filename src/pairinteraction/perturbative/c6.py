# SPDX-FileCopyrightText: 2025 Pairinteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

from typing import TYPE_CHECKING, Optional, Union, overload

from pairinteraction.perturbative.effective_system_pair import EffectiveSystemPair
from pairinteraction.units import QuantityScalar

if TYPE_CHECKING:
    from pairinteraction._wrapped.ket.ket_atom import KetAtom
    from pairinteraction.units import PintFloat


class C6(EffectiveSystemPair):
    def __init__(self, ket1: "KetAtom", ket2: "KetAtom") -> None:
        if ket1.species == ket2.species and ket1 != ket2:
            raise ValueError(
                "If you want to calculate 2nd order perturbations of two different states a and b "
                "(of the same species), "
                "please use the EffectiveSystemPair([(a,b), (b,a)]) class."
            )

        super().__init__([(ket1, ket2)])
        self.set_perturbation_order(2)

        # Set some default distance. The exact value does not matter for the C6 value,
        # but too small distances result in warnings about the basisvec admixtures)
        self.set_distance(100, unit="micrometer")

    @overload
    def get(self, unit: None = None) -> "PintFloat": ...

    @overload
    def get(self, unit: str) -> float: ...

    def get(self, unit: Optional[str] = None) -> Union[float, "PintFloat"]:
        """Get the C6 coefficient of the pair interaction between the specified ket1 and ket2.

        Args:
            unit: The unit in which to return the C6 coefficient.
                If None, returns a pint object.

        """
        h_eff_pint = self.get_effective_hamiltonian(return_order=2)
        distance = self.system_pair.get_distance()
        c6_pint = h_eff_pint[0, 0] * distance**6  # type: ignore [index]  # pint does not know it can be indexed
        return QuantityScalar.convert_pint_to_user(c6_pint, "c6", unit)
