# SPDX-FileCopyrightText: 2024 PairInteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

from __future__ import annotations

import threading
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from .utils import PairinteractionModule


def test_diagonalize_releases_gil(pi_module: PairinteractionModule) -> None:
    """Test that the main thread can proceed while another thread diagonalizes."""
    basis = pi_module.BasisAtom("Rb", n=(55, 68), l=(0, 4))
    systems = [
        pi_module.SystemAtom(basis).set_electric_field([0, 0, electric_field], unit="V/cm")
        for electric_field in np.linspace(0, 20, 25)
    ]

    def diagonalize_in_background() -> None:
        pi_module.diagonalize(systems, diagonalizer="eigen", sort_by_energy=True)

    worker = threading.Thread(target=diagonalize_in_background)
    worker.start()

    main_thread_progress = 0
    while worker.is_alive():
        main_thread_progress += 1

    worker.join()

    assert main_thread_progress > 0
