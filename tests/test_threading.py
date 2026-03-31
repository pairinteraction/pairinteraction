# SPDX-FileCopyrightText: 2024 PairInteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

from __future__ import annotations

import ctypes
import ctypes.util
import sys
import threading
import time
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from .utils import PairinteractionModule


def test_diagonalize_releases_gil(pi_module: PairinteractionModule) -> None:
    """Test that the main thread can proceed while another thread diagonalizes."""
    # Main thread processes during diagonalization because the GIL is released
    basis = pi_module.BasisAtom("Rb", n=(55, 68), l=(0, 2))
    systems = [
        pi_module.SystemAtom(basis).set_electric_field([0, 0, electric_field], unit="V/cm")
        for electric_field in np.linspace(0, 20, 10)
    ]

    barrier = threading.Barrier(2)

    def diagonalize_in_background() -> None:
        barrier.wait()
        pi_module.diagonalize(systems, diagonalizer="eigen", sort_by_energy=True)

    main_thread_progress_diagonalization = 0
    worker = threading.Thread(target=diagonalize_in_background)

    worker.start()
    barrier.wait()
    t0 = time.perf_counter()
    while worker.is_alive():
        main_thread_progress_diagonalization += 1
    worker.join()
    duration_diagonalization = time.perf_counter() - t0

    print(
        f"Time taken for diagonalization: {duration_diagonalization:.2f} seconds, "
        f"main thread progress: {main_thread_progress_diagonalization}"
    )

    # Main thread does not make significant progress during a GIL-holding C function
    if sys.platform == "win32":
        if not (kernel32_name := ctypes.util.find_library("kernel32")):
            raise RuntimeError("Could not find kernel32 library.")
        gil_holding_sleep = ctypes.PyDLL(kernel32_name).Sleep
        gil_holding_sleep.argtypes = [ctypes.c_uint32]
        sleep_duration = 250
    else:
        if not (libc_name := ctypes.util.find_library("c") or "libc.so.6"):
            raise RuntimeError("Could not find libc library.")
        gil_holding_sleep = ctypes.PyDLL(libc_name).usleep
        gil_holding_sleep.argtypes = [ctypes.c_uint]
        sleep_duration = 250_000

    barrier.reset()

    def c_function_holding_gil() -> None:
        barrier.wait()
        gil_holding_sleep(sleep_duration)

    main_thread_progress_python = 0
    worker = threading.Thread(target=c_function_holding_gil)

    worker.start()
    barrier.wait()
    t0 = time.perf_counter()
    while worker.is_alive():
        main_thread_progress_python += 1
    worker.join()
    duration_python = time.perf_counter() - t0

    print(
        f"Time taken for GIL-holding C function: {duration_python:.2f} seconds, "
        f"main thread progress: {main_thread_progress_python}"
    )

    assert (
        main_thread_progress_diagonalization / duration_diagonalization
        > 5 * main_thread_progress_python / duration_python
    ), (
        "The main thread did not make significant progress during diagonalization, "
        "suggesting that the GIL was not released."
    )
