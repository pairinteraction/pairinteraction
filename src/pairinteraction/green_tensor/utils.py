# SPDX-FileCopyrightText: 2024 PairInteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

# ruff: noqa: N802, N806

from __future__ import annotations

import math
from typing import TYPE_CHECKING, Callable, Literal, TypeVar

import numpy as np
import scipy.constants as const
from numba import njit
from scipy.integrate import quad

from pairinteraction.green_tensor.bessel_function import bessel_function_0, bessel_function_1, bessel_function_2

if TYPE_CHECKING:
    from pairinteraction.units import NDArray

    Entries = Literal["xx", "xy", "xz", "yx", "yy", "yz", "zx", "zy", "zz"]

    from typing_extensions import ParamSpec

    P = ParamSpec("P")
    R = TypeVar("R")

    def njit(cache: bool) -> Callable[[Callable[P, R]], Callable[P, R]]: ...


def green_tensor_homogeneous(
    r_a: NDArray, r_b: NDArray, omega: float, epsilon0: complex, *, only_real_part: bool = False
) -> NDArray:
    """Homogeneous Green Tensor for two atoms in cartesian coordinates in an infinite homogeneous medium.

    The function used is from equation 2 of the paper:
    "Dispersionless subradiant photon storage in one-dimensional emitter chains"

    Args:
        r_a: Position vector of atom A in meters
        r_b: Position vector of atom B in meters
        omega: Angular frequency in Hz
        epsilon0: Electric permittivity of the medium (dimensionless, complex)
        only_real_part: If True, only the real part of the Green tensor is calculated (default: False)

    Returns: The 3x3 Homogeneous Green Tensor (general complex values) (1/m)

    """
    k_vac = omega / const.c  # magnitude of wave vector in vacuum
    k0 = k_vac * np.sqrt(epsilon0)  # magnitude of wave vector in medium with permittivity epsilon0
    r = r_a - r_b
    r_norm = np.linalg.norm(r)

    prefactor = np.exp(1j * k0 * r_norm) / (4 * np.pi * k0**2 * r_norm**3)
    result: NDArray = prefactor * (
        (k0**2 * r_norm**2 + 1j * k0 * r_norm - 1) * np.eye(3)
        + (-(k0**2) * r_norm**2 - 3j * k0 * r_norm + 3) * np.outer(r, r) / r_norm**2
    )

    if only_real_part:
        return np.real(result)
    return result


@njit(cache=True)
def branch(epsilon: complex, k: float, k_rho: complex) -> complex:
    """Calculate the perpendicular wave vector component with positive imaginary part.

    Args:
        epsilon: Electric permittivity of the medium (dimensionless, complex)
        k: Magnitude of wave number in vacuum (1/m)
        k_rho: In-plane wave vector component (1/m)

    Returns: The perpendicular wave vector component (1/m)

    """
    return np.sqrt(epsilon * k**2 - k_rho**2 + 0j)


"""The following functions are used from Appendix B of the paper:
"Modified dipole-dipole interaction and dissipation in an atomic ensemble near surfaces"
and are needed to calculate the scattering Green Tensor for two atoms between two planar surfaces.
"""


@njit(cache=True)
def rs(kz: complex, k1z: complex) -> complex:
    """Calculate the Fresnel reflection coefficient for s-polarized light.

    Args:
        kz: Perpendicular wave vector component in the medium between the two surfaces (1/m)
        k1z: Perpendicular wave vector component in the upper or lower medium (1/m)

    Returns: The Fresnel reflection coefficient for s-polarized light (dimensionless, complex)

    """
    return (kz - k1z) / (kz + k1z)


@njit(cache=True)
def rp(kz: complex, k1z: complex, epsilon: complex) -> complex:
    """Calculate the Fresnel reflection coefficient for p-polarized light.

    Args:
        kz: Perpendicular wave vector component in the medium between the two surfaces (1/m)
        k1z: Perpendicular wave vector component in the upper or lower medium (1/m)
        epsilon: Electric permittivity of the upper or lower medium (dimensionless, complex)

    Returns: The Fresnel reflection coefficient for p-polarized light (dimensionless, complex)

    """
    return (epsilon * kz - k1z) / (epsilon * kz + k1z)


@njit(cache=True)
def D(r_plus: complex, r_minus: complex, kz: complex, h: float) -> complex:
    """Calculate the denominator term D used in the scattering Green Tensor matrix elements.

    Args:
        r_plus: Fresnel reflection coefficient for the upper surface (dimensionless, complex)
        r_minus: Fresnel reflection coefficient for the lower surface (dimensionless, complex)
        kz: Perpendicular wave vector component in the medium between the two surfaces (1/m)
        h: Distance between the two surfaces (m)

    Returns: The value of the denominator term D (dimensionless, complex)

    """
    return 1 - r_plus * r_minus * np.exp(2j * kz * h)


@njit(cache=True)
def A_plus(r_plus: complex, r_minus: complex, kz: complex, h: float, z_ges: float, z_ab: float) -> complex:
    """Calculate the numerator term A_plus used in the scattering Green Tensor matrix elements.

    Args:
        r_plus: Fresnel reflection coefficient for the upper surface (dimensionless, complex)
        r_minus: Fresnel reflection coefficient for the lower surface (dimensionless, complex)
        kz: Perpendicular wave vector component in the medium between the two surfaces (1/m)
        h: Distance between the two surfaces (m)
        z_ges: Total z-coordinate (m)
        z_ab: Difference in z-coordinates (m)

    Returns: The value of the numerator term A_plus (dimensionless, complex)

    """
    return (
        r_minus * np.exp(1j * kz * (z_ges - h))
        + r_plus * np.exp(-1j * kz * (z_ges - h))
        + 2 * r_plus * r_minus * np.cos(kz * z_ab) * np.exp(1j * kz * h)
    ) / D(r_plus, r_minus, kz, h)


@njit(cache=True)
def A_minus(r_plus: complex, r_minus: complex, kz: complex, h: float, z_ges: float, z_ab: float) -> complex:
    """Calculate the numerator term A_minus used in the scattering Green Tensor matrix elements.

    Args:
        r_plus: Fresnel reflection coefficient for the upper surface (dimensionless, complex)
        r_minus: Fresnel reflection coefficient for the lower surface (dimensionless, complex)
        kz: Perpendicular wave vector component in the medium between the two surfaces (1/m)
        h: Distance between the two surfaces (m)
        z_ges: Total z-coordinate (m)
        z_ab: Difference in z-coordinates (m)

    Returns: The value of the numerator term A_minus (dimensionless, complex)

    """
    return (
        r_minus * np.exp(1j * kz * (z_ges - h))
        + r_plus * np.exp(-1j * kz * (z_ges - h))
        - 2 * r_plus * r_minus * np.cos(kz * z_ab) * np.exp(1j * kz * h)
    ) / D(r_plus, r_minus, kz, h)


@njit(cache=True)
def B_plus(r_plus: complex, r_minus: complex, kz: complex, h: float, z_ges: float, z_ab: float) -> complex:
    """Calculate the numerator term B_plus used in the scattering Green Tensor matrix elements.

    Args:
        r_plus: Fresnel reflection coefficient for the upper surface (dimensionless, complex)
        r_minus: Fresnel reflection coefficient for the lower surface (dimensionless, complex)
        kz: Perpendicular wave vector component in the medium between the two surfaces (1/m)
        h: Distance between the two surfaces (m)
        z_ges: Total z-coordinate (m)
        z_ab: Difference in z-coordinates (m)

    Returns: The value of the numerator term B_plus (dimensionless, complex)

    """
    return (
        r_minus * np.exp(1j * kz * (z_ges - h))
        + r_plus * np.exp(-1j * kz * (z_ges - h))
        + 2j * r_plus * r_minus * np.sin(kz * z_ab) * np.exp(1j * kz * h)
    ) / D(r_plus, r_minus, kz, h)


@njit(cache=True)
def B_minus(r_plus: complex, r_minus: complex, kz: complex, h: float, z_ges: float, z_ab: float) -> complex:
    """Calculate the numerator term B_minus used in the scattering Green Tensor matrix elements.

    Args:
        r_plus: Fresnel reflection coefficient for the upper surface (dimensionless, complex)
        r_minus: Fresnel reflection coefficient for the lower surface (dimensionless, complex)
        kz: Perpendicular wave vector component in the medium between the two surfaces (1/m)
        h: Distance between the two surfaces (m)
        z_ges: Total z-coordinate (m)
        z_ab: Difference in z-coordinates (m)

    Returns: The value of the numerator term B_minus (dimensionless, complex)

    """
    return (
        r_minus * np.exp(1j * kz * (z_ges - h))
        + r_plus * np.exp(-1j * kz * (z_ges - h))
        - 2j * r_plus * r_minus * np.sin(kz * z_ab) * np.exp(1j * kz * h)
    ) / D(r_plus, r_minus, kz, h)


def Gs(
    kz: complex,
    h: float,
    k_rho: complex,
    rho: float,
    phi: float,
    rs_plus: complex,
    rs_minus: complex,
    z_ges: float,
    z_ab: float,
    entry: Entries,
) -> complex:
    """Calculate the Gs part of the scattering Green Tensor."""
    if entry in ["xz", "yz", "zx", "zy", "zz"]:
        return 0
    As_plus = A_plus(rs_plus, rs_minus, kz, h, z_ges, z_ab)
    J2 = bessel_function_2(k_rho * rho)
    if entry in ["xy", "yx"]:
        return -As_plus / 2 * J2 * math.sin(2 * phi)
    J0 = bessel_function_0(k_rho * rho)
    if entry == "xx":
        return As_plus / 2 * (J0 + J2 * math.cos(2 * phi))
    if entry == "yy":
        return As_plus / 2 * (J0 - J2 * math.cos(2 * phi))

    raise ValueError(f"Invalid entry '{entry}' for Gs function.")


def Gp(  # noqa: PLR0911
    kz: complex,
    h: float,
    k_rho: complex,
    rho: float,
    phi: float,
    rp_plus: complex,
    rp_minus: complex,
    z_ges: float,
    z_ab: float,
    entry: Entries,
) -> complex:
    """Calculate the Gp part of the scattering Green Tensor."""
    if entry in ["xz", "zx", "yz", "zy"]:
        J1 = bessel_function_1(k_rho * rho)
        if entry == "zx":
            Bp_minus = B_minus(rp_plus, rp_minus, kz, h, z_ges, z_ab)
            return -1j * (k_rho / kz) * Bp_minus * J1 * math.cos(phi)
        Bp_plus = B_plus(rp_plus, rp_minus, kz, h, z_ges, z_ab)
        if entry == "xz":
            return 1j * (k_rho / kz) * Bp_plus * J1 * math.cos(phi)
        if entry == "yz":
            return 1j * (k_rho / kz) * Bp_plus * J1 * math.sin(phi)
        if entry == "zy":
            return -1j * (k_rho / kz) * Bp_plus * J1 * math.sin(phi)

    if entry == "zz":
        J0 = bessel_function_0(k_rho * rho)
        Ap_plus = A_plus(rp_plus, rp_minus, kz, h, z_ges, z_ab)
        return -(k_rho**2 / kz**2) * Ap_plus * J0

    J2 = bessel_function_2(k_rho * rho)
    Ap_minus = A_minus(rp_plus, rp_minus, kz, h, z_ges, z_ab)
    if entry in ["xy", "yx"]:
        return Ap_minus / 2 * J2 * math.sin(2 * phi)

    J0 = bessel_function_0(k_rho * rho)
    if entry == "xx":
        return Ap_minus / 2 * (J0 - J2 * math.cos(2 * phi))
    if entry == "yy":
        return Ap_minus / 2 * (J0 + J2 * math.cos(2 * phi))

    raise ValueError(f"Invalid entry '{entry}' for Gp function.")


""" The integrals for the scattering Green Tensor are evaluated in two parts:
    The first part is along an elliptical path from 0 to 2*k_maj in the complex plane to avoid singularities,
    and the second part is along the real axis from 2*k_maj to an upper limit.

    The methods for evaluating these integrals are explained in the papers:
    - "Accurate and efficient computation of the Green's tensor for stratified media" section III.A
    - "Challenges in Computational Electromagnetics: Analysis and Optimization of Planar Multilayered Structures"
        section 2.4.1
    """


@njit(cache=True)
def integrand_ellipse_partial(
    t: complex,
    k_maj: float,
    k_min: float,
    k0: float,
    epsilon1: complex,
    epsilon2: complex,
    h: float,
) -> tuple[complex, complex, complex, complex, complex, complex, complex]:
    # elliptical path, substitution
    k_rho = k_maj * (1 + np.cos(t)) - 1j * k_min * np.sin(t)
    dk_rho = -k_maj * np.sin(t) - 1j * k_min * np.cos(t)

    # Wave vector components
    kz = branch(1, k0, k_rho)
    k1z = branch(epsilon1, k0, k_rho)
    k2z = branch(epsilon2, k0, k_rho)

    rs_plus = rs(kz, k1z)
    rs_minus = rs(kz, k2z)
    rp_plus = rp(kz, k1z, epsilon1)
    rp_minus = rp(kz, k2z, epsilon2)

    prefactor = 1j / (4 * np.pi) * (k_rho / kz) * np.exp(1j * kz * h) * dk_rho
    return k_rho, kz, rs_plus, rs_minus, rp_plus, rp_minus, prefactor


def integrand_ellipse(
    t: complex,
    k_maj: float,
    k_min: float,
    k0: float,
    epsilon1: complex,
    epsilon2: complex,
    h: float,
    rho: float,
    phi: float,
    z_ges: float,
    z_ab: float,
    entry: Entries,
    real_or_imag: str,
) -> complex:
    k_rho, kz, rs_plus, rs_minus, rp_plus, rp_minus, prefactor = integrand_ellipse_partial(
        t, k_maj, k_min, k0, epsilon1, epsilon2, h
    )

    gs = Gs(kz, h, k_rho, rho, phi, rs_plus, rs_minus, z_ges, z_ab, entry)
    gp = Gp(kz, h, k_rho, rho, phi, rp_plus, rp_minus, z_ges, z_ab, entry)

    integrand = prefactor * (gs - (kz**2 / k0**2) * gp)
    if real_or_imag == "real":
        return np.real(integrand)
    if real_or_imag == "imag":
        return np.imag(integrand)
    raise ValueError("real_or_imag must be 'real' or 'imag'")


def elliptic_integral(
    omega: float,
    h: float,
    rho: float,
    phi: float,
    epsilon0: complex,
    epsilon1: complex,
    epsilon2: complex,
    z_ges: float,
    z_ab: float,
    entry: Entries,
    *,
    only_real_part: bool = False,
) -> complex:
    """Evaluate the elliptic part of the integral.

    Args:
        omega: Angular frequency in Hz
        h: Distance between the two surfaces in meters
        rho: In-plane distance between the two atoms in meters
        phi: Angle between the in-plane distance vector and the x-axis in radians
        epsilon0: Electric permittivity of the medium between the two surfaces (dimensionless, complex)
        epsilon1: Electric permittivity of the upper medium (dimensionless, complex)
        epsilon2: Electric permittivity of the lower medium (dimensionless, complex)
        z_ges: Sum of the z-positions of the two atoms in meters
        z_ab: Difference of the z-positions of the two atoms in meters
        entry: Entry of the Green tensor to calculate
        only_real_part: If True, only the real part of the integral is calculated (default: False)

    Returns: The value of the integral along the elliptical path as a complex number (1/m)

    """
    k_vac = omega / const.c  # magnitude of wave vector in vacuum
    k0 = k_vac * np.sqrt(epsilon0)

    # Elliptical path in complex plane to avoid singularities (Integral from 0 to 2k_maj)
    k1 = k_vac * np.sqrt(epsilon1)
    k2 = k_vac * np.sqrt(epsilon2)
    kl_max = max(np.real(k0), np.real(k1), np.real(k2))

    k_maj = (kl_max + k_vac) / 2  # major axis of ellipse
    k_min = min(k_vac, 1 / rho) if rho != 0 else k_vac

    args = (k_maj, k_min, k0, epsilon1, epsilon2, h, rho, phi, z_ges, z_ab, entry)

    real_ellipse, _ = quad(integrand_ellipse, np.pi, 0, args=(*args, "real"), epsrel=1e-9, limit=1000)  # type: ignore [arg-type]
    if only_real_part:
        return real_ellipse
    imag_ellipse, _ = quad(integrand_ellipse, np.pi, 0, args=(*args, "imag"), epsrel=1e-9, limit=1000)  # type: ignore [arg-type]
    return real_ellipse + 1j * imag_ellipse


def integrand_real(
    k_rho: complex,
    k0: float,
    epsilon1: complex,
    epsilon2: complex,
    h: float,
    rho: float,
    phi: float,
    z_ges: float,
    z_ab: float,
    entry: Entries,
    real_or_imag: str,
) -> complex:
    kz = branch(1, k0, k_rho)
    k1z = branch(epsilon1, k0, k_rho)
    k2z = branch(epsilon2, k0, k_rho)

    rs_plus = rs(kz, k1z)
    rs_minus = rs(kz, k2z)
    rp_plus = rp(kz, k1z, epsilon1)
    rp_minus = rp(kz, k2z, epsilon2)

    # Integrand
    integrand = (
        1j
        / (4 * np.pi)
        * (
            Gs(kz, h, k_rho, rho, phi, rs_plus, rs_minus, z_ges, z_ab, entry)
            - (kz**2 / k0**2) * Gp(kz, h, k_rho, rho, phi, rp_plus, rp_minus, z_ges, z_ab, entry)
        )
        * (k_rho / kz)
        * np.exp(1j * kz * h)
    )
    if real_or_imag == "real":
        return np.real(integrand)
    if real_or_imag == "imag":
        return np.imag(integrand)
    raise ValueError("real_or_imag must be 'real' or 'imag'")


def real_axis_integral(
    omega: float,
    h: float,
    rho: float,
    phi: float,
    epsilon0: complex,
    epsilon1: complex,
    epsilon2: complex,
    z_ges: float,
    z_ab: float,
    entry: Entries,
    upper_limit: float,
    *,
    only_real_part: bool = False,
) -> complex:
    """Evaluate the real axis part of the integral.

    Args:
        omega: Angular frequency in Hz
        h: Distance between the two surfaces in meters
        rho: In-plane distance between the two atoms in meters
        phi: Angle between the in-plane distance vector and the x-axis in radians
        epsilon0: Electric permittivity of the medium between the two surfaces (dimensionless, complex)
        epsilon1: Electric permittivity of the upper medium (dimensionless, complex)
        epsilon2: Electric permittivity of the lower medium (dimensionless, complex)
        z_ges: Sum of the z-positions of the two atoms in meters
        z_ab: Difference of the z-positions of the two atoms in meters
        entry: Entry of the Green tensor to calculate
        upper_limit: Upper limit for the real axis integral (1/m)
        only_real_part: If True, only the real part of the integral is calculated (default: False)

    Returns: The value of the integral along the real axis as a complex number (1/m)

    """
    k_vac = omega / const.c  # magnitude of wave vector in vacuum
    k0 = k_vac * np.sqrt(epsilon0)
    k1 = k_vac * np.sqrt(epsilon1)
    k2 = k_vac * np.sqrt(epsilon2)

    kl_max = max(np.real(k0), np.real(k1), np.real(k2))
    k_maj = (kl_max + k_vac) / 2

    args = (k0, epsilon1, epsilon2, h, rho, phi, z_ges, z_ab, entry)

    real_real, _ = quad(
        integrand_real,  # type: ignore [arg-type]
        2 * k_maj,
        upper_limit * np.real(k0),
        args=(*args, "real"),  # type: ignore [arg-type]
        limit=1000,
        epsrel=1e-9,
    )
    if only_real_part:
        return real_real
    imag_real, _ = quad(
        integrand_real,  # type: ignore [arg-type]
        2 * k_maj,
        upper_limit * np.real(k0),
        args=(*args, "imag"),  # type: ignore [arg-type]
        limit=1000,
        epsrel=1e-9,
    )
    return real_real + 1j * imag_real


def green_tensor_scattered(
    r_a: NDArray,
    r_b: NDArray,
    omega: float,
    epsilon0: complex,
    epsilon1: complex,
    epsilon2: complex,
    h: float,
    *,
    only_real_part: bool = False,
) -> NDArray:
    """Assemble the total scattering Green tensor.

    Args:
        r_a: Position vector of atom A (m)
        r_b: Position vector of atom B (m)
        omega: Angular frequency in (1/s)
        epsilon0: Electric permittivity of the medium between the two surfaces (dimensionless, complex)
        epsilon1: Electric permittivity of the upper medium (dimensionless, complex)
        epsilon2: Electric permittivity of the lower medium (dimensionless, complex)
        h: Distance between the two surfaces (m)
        only_real_part: If True, only the real part of the Green tensor is calculated (default: False)

    Returns: The 3x3 Scattering Green Tensor (general complex values) (1/m)

    """
    r = r_a - r_b
    rho = np.sqrt(r[0] ** 2 + r[1] ** 2)
    z_alpha = r_a[2]
    z_beta = r_b[2]
    z_ges = z_alpha + z_beta
    z_ab = r[2]
    phi = np.arccos(r[0] / rho) if rho != 0 else 0

    # Estimate the upper limit for the real axis integral
    k_vac = omega / const.c  # magnitude of wave vector in vacuum
    k0 = k_vac * np.sqrt(epsilon0)
    upper_limit_real_integral = np.sqrt((745 / (np.real(k0) * h)) ** 2 + 1)

    gt_total = np.zeros((3, 3), dtype=complex)
    for i, ix in enumerate(["x", "y", "z"]):
        for j, jx in enumerate(["x", "y", "z"]):
            entry: Entries = ix + jx  # type: ignore [assignment]
            g_ij_elliptic = elliptic_integral(
                omega, h, rho, phi, epsilon0, epsilon1, epsilon2, z_ges, z_ab, entry, only_real_part=only_real_part
            )
            g_ij_real = real_axis_integral(
                omega,
                h,
                rho,
                phi,
                epsilon0,
                epsilon1,
                epsilon2,
                z_ges,
                z_ab,
                entry,
                upper_limit_real_integral,
                only_real_part=only_real_part,
            )
            gt_total[i][j] = g_ij_elliptic + g_ij_real

    return gt_total


def green_tensor_total(
    r_a: NDArray,
    r_b: NDArray,
    omega: float,
    epsilon0: complex,
    epsilon1: complex,
    epsilon2: complex,
    h: float,
    *,
    only_real_part: bool = False,
) -> NDArray:
    """Assemble the total Green tensor.

    Args:
        r_a: Position vector of atom A (m)
        r_b: Position vector of atom B (m)
        omega: Angular frequency in (1/s)
        epsilon0: Electric permittivity of the medium between the two surfaces (dimensionless, complex)
        epsilon1: Electric permittivity of the upper medium (dimensionless, complex)
        epsilon2: Electric permittivity of the lower medium (dimensionless, complex)
        h: Distance between the two surfaces (m).
        only_real_part: If True, only the real part of the Green tensor is calculated (default: False)

    Returns: The 3x3 Total Green Tensor (general complex values) (1/m)

    """
    gt_scat = green_tensor_scattered(r_a, r_b, omega, epsilon0, epsilon1, epsilon2, h, only_real_part=only_real_part)
    gt_hom = green_tensor_homogeneous(r_a, r_b, omega, epsilon0, only_real_part=only_real_part)

    return gt_scat + gt_hom
