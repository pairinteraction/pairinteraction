# SPDX-FileCopyrightText: 2025 Pairinteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

import logging
from typing import TYPE_CHECKING, Optional, Union, overload

from pairinteraction.perturbative.effective_hamiltonian import EffectiveHamiltonian
from pairinteraction.units import QuantityScalar

if TYPE_CHECKING:
    from pairinteraction._wrapped.ket.ket_atom import KetAtom
    from pairinteraction.units import PintFloat

logger = logging.getLogger(__name__)


class C3(EffectiveHamiltonian):
    def __init__(self, ket1: "KetAtom", ket2: "KetAtom") -> None:
        super().__init__([(ket1, ket2), (ket2, ket1)])
        self.set_perturbation_order(1)
        self.set_distance(20, unit="micrometer")

    @overload
    def get_c3(self, unit: None = None) -> "PintFloat": ...

    @overload
    def get_c3(self, unit: str) -> float: ...

    def get_c3(self, unit: Optional[str] = None) -> Union[float, "PintFloat"]:
        """Get the C3 coefficient of the pair interaction between the specified ket1 and ket2.

        Args:
            unit: The unit in which to return the C3 coefficient.
                If None, returns a pint object.

        """
        h_eff_pint = self.get_effective_hamiltonian(return_order=1)
        distance = self.system_pair.get_distance()
        c3_pint = h_eff_pint[1, 0] * distance**3  # type: ignore [index]  # pint does not know it can be indexed
        return QuantityScalar.convert_pint_to_user(c3_pint, "c3", unit)
