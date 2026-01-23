# SPDX-FileCopyrightText: 2025 PairInteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

from functools import lru_cache

from scipy.special import jv as scipy_bessel_function


@lru_cache(maxsize=100_000)
def bessel_function_0(z: complex) -> complex:
    return scipy_bessel_function(0, z)


@lru_cache(maxsize=100_000)
def bessel_function_1(z: complex) -> complex:
    return scipy_bessel_function(1, z)


@lru_cache(maxsize=100_000)
def bessel_function_2(z: complex) -> complex:
    return scipy_bessel_function(2, z)
    # for even a bit more speedup:
    # return (2.0 / z) * bessel_function_1(z) - bessel_function_0(z)  # noqa: ERA001
