# SPDX-FileCopyrightText: 2024 PairInteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

import numpy as np
import scipy.constants as const
from numpy.typing import NDArray
from scipy.integrate import quad
from scipy.special import jn

# TODO implement utils for green tensors

def green_tensor_homogeneous(r_a: NDArray, r_b: NDArray, omega: float, epsilon0: complex) -> NDArray:
    """ Homogeneous Green Tensor for two atoms in cartesian coordinates in an infinite homogeneous medium.

    The function used is from equation 2 of the paper: "Dispersionless subradiant photon storage in one-dimensional emitter chains"

    Args:
        r_a: Position vector of atom A in meters
        r_b: Position vector of atom B in meters
        omega: Angular frequency in 1/seconds
        epsilon0: Electric permittivity of the medium (dimensionless, complex)

    Returns: The 3x3 Homogeneous Green Tensor (general complex values) (1/m)
    """
    k_vac = omega / const.c  # magnitude of wave vector in vacuum
    k0 = k_vac * np.sqrt(epsilon0 + 0j)  # magnitude of wave vector in medium with permittivity epsilon0
    r = r_a - r_b
    r_norm = np.linalg.norm(r)

    prefactor = np.exp(1j * k0 * r_norm) / (4 * np.pi * k0**2 * r_norm**3)
    gt_hom = prefactor * ((k0**2 * r_norm**2 + 1j * k0 * r_norm - 1) * np.eye(3) + (-(k0**2) * r_norm**2 - 3j * k0 * r_norm + 3) * np.outer(r, r) / r_norm**2)

    return gt_hom


def branch(epsilon, k, k_rho):
    """ Calculate the perpendicular wave vector component with positive imaginary part.

    Args:
        epsilon: Electric permittivity of the medium (dimensionless, complex)
        k: Magnitude of wave number in vacuum (1/m)
        k_rho: In-plane wave vector component (1/m)

    Returns: The perpendicular wave vector component (1/m)
    """
    k_perp = np.sqrt(epsilon * k**2 - k_rho**2 + 0j)
    return k_perp


""" The following functions are used from Appendix B of the paper: "Modified dipole-dipole interaction and dissipation in an atomic ensemble near surfaces"
    and are needed to calculate the scattering Green Tensor for two atoms between two planar surfaces."""


def rs(kz: complex, k1z: complex) -> complex:
    """ Calculate the Fresnel reflection coefficient for s-polarized light.

    Args:
        kz: Perpendicular wave vector component in the medium between the two surfaces (1/m)
        k1z: Perpendicular wave vector component in the upper or lower medium (1/m)

    Returns: The Fresnel reflection coefficient for s-polarized light (dimensionless, complex)
    """
    return (kz - k1z) / (kz + k1z)


def rp(kz: complex, k1z: complex, epsilon: complex) -> complex:
    """ Calculate the Fresnel reflection coefficient for p-polarized light.

    Args:
        kz: Perpendicular wave vector component in the medium between the two surfaces (1/m)
        k1z: Perpendicular wave vector component in the upper or lower medium (1/m)
        epsilon: Electric permittivity of the upper or lower medium (dimensionless, complex)

    Returns: The Fresnel reflection coefficient for p-polarized light (dimensionless, complex)
    """
    return (epsilon * kz - k1z) / (epsilon * kz + k1z)


# Prefactors for matrix elements of scattering Green's Tensor
def D(r_plus: complex, r_minus: complex, kz: complex, h: float) -> complex:
    """ Calculate the denominator term D used in the scattering Green Tensor matrix elements.

    Args:
        r_plus: Fresnel reflection coefficient for the upper surface (dimensionless, complex)
        r_minus: Fresnel reflection coefficient for the lower surface (dimensionless, complex)
        kz: Perpendicular wave vector component in the medium between the two surfaces (1/m)
        h: Distance between the two surfaces (m)

    Returns: The value of the denominator term D (dimensionless, complex)
    """
    return 1 - r_plus * r_minus * np.exp(2j * kz * h)


def A_plus(r_plus: complex, r_minus: complex, kz: complex, h: float, z_ges: float, z_ab: float) -> complex:
    """ Calculate the numerator term A_plus used in the scattering Green Tensor matrix elements.

    Args:
        r_plus: Fresnel reflection coefficient for the upper surface (dimensionless, complex)
        r_minus: Fresnel reflection coefficient for the lower surface (dimensionless, complex)
        kz: Perpendicular wave vector component in the medium between the two surfaces (1/m)
        h: Distance between the two surfaces (m)
        z_ges: Total z-coordinate (m)
        z_AB: Difference in z-coordinates (m)

    Returns: The value of the numerator term A_plus (dimensionless, complex)
    """
    return (
        r_minus * np.exp(1j * kz * (z_ges - h))
        + r_plus * np.exp(-1j * kz * (z_ges - h))
        + 2 * r_plus * r_minus * np.cos(kz * z_ab) * np.exp(1j * kz * h)
    ) / D(r_plus, r_minus, kz, h)


def A_minus(r_plus: complex, r_minus: complex, kz: complex, h: float, z_ges: float, z_ab: float) -> complex:
    """ Calculate the numerator term A_minus used in the scattering Green Tensor matrix elements.

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


def B_plus(r_plus: complex, r_minus: complex, kz: complex, h: float, z_ges: float, z_ab: float) -> complex:
    """ Calculate the numerator term B_plus used in the scattering Green Tensor matrix elements.

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


def B_minus(r_plus: complex, r_minus: complex, kz: complex, h: float, z_ges: float, z_ab: float) -> complex:
    """ Calculate the numerator term B_minus used in the scattering Green Tensor matrix elements.

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


''' Calculation of matrix elements of the scattering Green's Tensor '''
def Gs_xx(
    kz: complex,
    h: float,
    k_rho: complex,
    rho: float,
    phi: float,
    rs_plus: complex,
    rs_minus: complex,
    rp_plus: complex,
    rp_minus: complex,
    z_ges: float,
    z_ab: float,
) -> complex:
    # Bessel functions
    J0 = jn(0, k_rho * rho)
    J2 = jn(2, k_rho * rho)

    As_plus = A_plus(rs_plus, rs_minus, kz, h, z_ges, z_ab)

    return As_plus / 2 * (J0 + J2 * np.cos(2 * phi))


def Gs_yy(
    kz: complex,
    h: float,
    k_rho: complex,
    rho: float,
    phi: float,
    rs_plus: complex,
    rs_minus: complex,
    rp_plus: complex,
    rp_minus: complex,
    z_ges: float,
    z_ab: float,
) -> complex:
    # Bessel functions
    J0 = jn(0, k_rho * rho)
    J2 = jn(2, k_rho * rho)

    As_plus = A_plus(rs_plus, rs_minus, kz, h, z_ges, z_ab)

    return As_plus / 2 * (J0 - J2 * np.cos(2 * phi))


def Gs_zz(
    kz: complex,
    h: float,
    k_rho: complex,
    rho: float,
    phi: float,
    rs_plus: complex,
    rs_minus: complex,
    rp_plus: complex,
    rp_minus: complex,
    z_ges: float,
    z_ab: float,
) -> complex:
    return 0


def Gp_xx(
    kz: complex,
    h: float,
    k_rho: complex,
    rho: float,
    phi: float,
    rs_plus: complex,
    rs_minus: complex,
    rp_plus: complex,
    rp_minus: complex,
    z_ges: float,
    z_ab: float,
) -> complex:
    # Bessel functions
    J0 = jn(0, k_rho * rho)
    J2 = jn(2, k_rho * rho)

    Ap_minus = A_minus(rp_plus, rp_minus, kz, h, z_ges, z_ab)

    return Ap_minus / 2 * (J0 - J2 * np.cos(2 * phi))


def Gp_yy(
    kz: complex,
    h: float,
    k_rho: complex,
    rho: float,
    phi: float,
    rs_plus: complex,
    rs_minus: complex,
    rp_plus: complex,
    rp_minus: complex,
    z_ges: float,
    z_ab: float,
) -> complex:
    # Bessel functions
    J0 = jn(0, k_rho * rho)
    J2 = jn(2, k_rho * rho)

    Ap_minus = A_minus(rp_plus, rp_minus, kz, h, z_ges, z_ab)

    return Ap_minus / 2 * (J0 + J2 * np.cos(2 * phi))


def Gp_zz(
    kz: complex,
    h: float,
    k_rho: complex,
    rho: float,
    phi: float,
    rs_plus: complex,
    rs_minus: complex,
    rp_plus: complex,
    rp_minus: complex,
    z_ges: float,
    z_ab: float,
) -> complex:
    # Bessel function
    J0 = jn(0, k_rho * rho)

    Ap_plus = A_plus(rp_plus, rp_minus, kz, h, z_ges, z_ab)

    return -(k_rho**2 / kz**2) * Ap_plus * J0


def Gs_xy(
    kz: complex,
    h: float,
    k_rho: complex,
    rho: float,
    phi: float,
    rs_plus: complex,
    rs_minus: complex,
    rp_plus: complex,
    rp_minus: complex,
    z_ges: float,
    z_ab: float,
) -> complex:
    J2 = jn(2, k_rho * rho)
    As_plus = A_plus(rs_plus, rs_minus, kz, h, z_ges, z_ab)

    return -As_plus / 2 * J2 * np.sin(2 * phi)


def Gs_yx(
    kz: complex,
    h: float,
    k_rho: complex,
    rho: float,
    phi: float,
    rs_plus: complex,
    rs_minus: complex,
    rp_plus: complex,
    rp_minus: complex,
    z_ges: float,
    z_ab: float,
) -> complex:
    J2 = jn(2, k_rho * rho)
    As_plus = A_plus(rs_plus, rs_minus, kz, h, z_ges, z_ab)

    return -As_plus / 2 * J2 * np.sin(2 * phi)


def Gs_xz(
    kz: complex,
    h: float,
    k_rho: complex,
    rho: float,
    phi: float,
    rs_plus: complex,
    rs_minus: complex,
    rp_plus: complex,
    rp_minus: complex,
    z_ges: float,
    z_ab: float,
) -> complex:
    return 0


def Gs_zx(
    kz: complex,
    h: float,
    k_rho: complex,
    rho: float,
    phi: float,
    rs_plus: complex,
    rs_minus: complex,
    rp_plus: complex,
    rp_minus: complex,
    z_ges: float,
    z_ab: float,
) -> complex:
    return 0


def Gs_yz(
    kz: complex,
    h: float,
    k_rho: complex,
    rho: float,
    phi: float,
    rs_plus: complex,
    rs_minus: complex,
    rp_plus: complex,
    rp_minus: complex,
    z_ges: float,
    z_ab: float,
) -> complex:
    return 0


def Gs_zy(
    kz: complex,
    h: float,
    k_rho: complex,
    rho: float,
    phi: float,
    rs_plus: complex,
    rs_minus: complex,
    rp_plus: complex,
    rp_minus: complex,
    z_ges: float,
    z_ab: float,
) -> complex:
    return 0


def Gp_xy(
    kz: complex,
    h: float,
    k_rho: complex,
    rho: float,
    phi: float,
    rs_plus: complex,
    rs_minus: complex,
    rp_plus: complex,
    rp_minus: complex,
    z_ges: float,
    z_ab: float,
) -> complex:
    J2 = jn(2, k_rho * rho)
    Ap_minus = A_minus(rp_plus, rp_minus, kz, h, z_ges, z_ab)

    return Ap_minus / 2 * J2 * np.sin(2 * phi)


def Gp_xz(
    kz: complex,
    h: float,
    k_rho: complex,
    rho: float,
    phi: float,
    rs_plus: complex,
    rs_minus: complex,
    rp_plus: complex,
    rp_minus: complex,
    z_ges: float,
    z_ab: float,
) -> complex:
    J1 = jn(1, k_rho * rho)
    Bp_plus = B_plus(rp_plus, rp_minus, kz, h, z_ges, z_ab)

    return 1j * (k_rho / kz) * Bp_plus * J1 * np.cos(phi)


def Gp_yx(
    kz: complex,
    h: float,
    k_rho: complex,
    rho: float,
    phi: float,
    rs_plus: complex,
    rs_minus: complex,
    rp_plus: complex,
    rp_minus: complex,
    z_ges: float,
    z_ab: float,
) -> complex:
    J2 = jn(2, k_rho * rho)
    Ap_minus = A_minus(rp_plus, rp_minus, kz, h, z_ges, z_ab)

    return Ap_minus / 2 * J2 * np.sin(2 * phi)


def Gp_yz(
    kz: complex,
    h: float,
    k_rho: complex,
    rho: float,
    phi: float,
    rs_plus: complex,
    rs_minus: complex,
    rp_plus: complex,
    rp_minus: complex,
    z_ges: float,
    z_ab: float,
) -> complex:
    J1 = jn(1, k_rho * rho)
    Bp_plus = B_plus(rp_plus, rp_minus, kz, h, z_ges, z_ab)

    return 1j * (k_rho / kz) * Bp_plus * J1 * np.sin(phi)


def Gp_zx(
    kz: complex,
    h: float,
    k_rho: complex,
    rho: float,
    phi: float,
    rs_plus: complex,
    rs_minus: complex,
    rp_plus: complex,
    rp_minus: complex,
    z_ges: float,
    z_ab: float,
) -> complex:
    J1 = jn(1, k_rho * rho)
    Bp_minus = B_minus(rp_plus, rp_minus, kz, h, z_ges, z_ab)

    return -1j * (k_rho / kz) * Bp_minus * J1 * np.cos(phi)


def Gp_zy(
    kz: complex,
    h: float,
    k_rho: complex,
    rho: float,
    phi: float,
    rs_plus: complex,
    rs_minus: complex,
    rp_plus: complex,
    rp_minus: complex,
    z_ges: float,
    z_ab: float,
) -> complex:
    J1 = jn(1, k_rho * rho)
    Bp_plus = B_plus(rp_plus, rp_minus, kz, h, z_ges, z_ab)

    return -1j * (k_rho / kz) * Bp_plus * J1 * np.sin(phi)

""" The integrals for the scattering Green Tensor are evaluated in two parts:
    The first part is along an elliptical path from 0 to 2*k_maj in the complex plane to avoid singularities,
    and the second part is along the real axis from 2*k_maj to an upper limit.

    The methods for evaluating these integrals are explained in the papers:
    "Accurate and efficient computation of the Green’s tensor for stratified media" section III.A
    "Challenges in Computational Electromagnetics: Analysis and Optimization of Planar Multilayered Structures" section 2.4.1
    """

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
    Gs,
    Gp,
) -> complex:
    """ Function to evaluate the elliptic part of the integral
    Args:
        omega: Angular frequency in 1/seconds
        h: Distance between the two surfaces in meters
        rho: In-plane distance between the two atoms in meters
        phi: Angle between the in-plane distance vector and the x-axis in radians
        epsilon0: Electric permittivity of the medium between the two surfaces (dimensionless, complex)
        epsilon1: Electric permittivity of the upper medium (dimensionless, complex)
        epsilon2: Electric permittivity of the lower medium (dimensionless, complex)
        z_ges: Sum of the z-positions of the two atoms in meters
        z_AB: Difference of the z-positions of the two atoms in meters
        gs: Function to calculate the s-polarized contribution to the Green tensor
        gp: Function to calculate the p-polarized contribution to the Green tensor

    Returns: The value of the integral along the elliptical path as a complex number (1/m)
    """
    k_vac = omega / const.c  # magnitude of wave vector in vacuum
    k0 = k_vac * np.sqrt(epsilon0 + 0j)

    # Elliptical path in complex plane to avoid singularities (Integral from 0 to 2k_maj)
    k1 = k_vac * np.sqrt(epsilon1 + 0j)
    k2 = k_vac * np.sqrt(epsilon2 + 0j)
    kl_max = max(np.real(k0), np.real(k1), np.real(k2))

    k_maj = (kl_max + k_vac) / 2  # major axis of ellipse
    k_min = (
        min(k_vac, 1 / rho) if rho != 0 else k_vac
    )

    def integrand_ellipse(t):
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

        # Integrand
        to_integrate = (
            1j
            / (4 * np.pi)
            * (
                Gs(kz, h, k_rho, rho, phi, rs_plus, rs_minus, rp_plus, rp_minus, z_ges, z_ab)
                - (kz**2 / k0**2) * Gp(kz, h, k_rho, rho, phi, rs_plus, rs_minus, rp_plus, rp_minus, z_ges, z_ab)
            )
            * (k_rho / kz)
            * np.exp(1j * kz * h)
            * dk_rho
        )

        return to_integrate

    # Real and imaginary part
    real_part_ellipse = lambda t: np.real(integrand_ellipse(t))
    imag_part_ellipse = lambda t: np.imag(integrand_ellipse(t))

    Re_ellipse, _ = quad(real_part_ellipse, np.pi, 0, epsrel=1e-9, limit=1000)
    Im_ellipse, _ = quad(imag_part_ellipse, np.pi, 0, epsrel=1e-9, limit=1000)

    return Re_ellipse + 1j * Im_ellipse

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
    Gs,
    Gp,
    upper_limit: float,
) -> complex:
    """ Function to evaluate the real axis part of the integral

    Args:
        omega: Angular frequency in 1/seconds
        h: Distance between the two surfaces in meters
        rho: In-plane distance between the two atoms in meters
        phi: Angle between the in-plane distance vector and the x-axis in radians
        epsilon0: Electric permittivity of the medium between the two surfaces (dimensionless, complex)
        epsilon1: Electric permittivity of the upper medium (dimensionless, complex)
        epsilon2: Electric permittivity of the lower medium (dimensionless, complex)
        z_ges: Sum of the z-positions of the two atoms in meters
        z_ab: Difference of the z-positions of the two atoms in meters
        gs: Function to calculate the s-polarized contribution to the Green tensor
        gp: Function to calculate the p-polarized contribution to the Green tensor
        upper_limit: Upper limit for the real axis integral (1/m)

    Returns: The value of the integral along the real axis as a complex number (1/m)
    """
    k_vac = omega / const.c  # magnitude of wave vector in vacuum
    k0 = k_vac * np.sqrt(epsilon0 + 0j)
    k1 = k_vac * np.sqrt(epsilon1 + 0j)
    k2 = k_vac * np.sqrt(epsilon2 + 0j)

    kl_max = max(np.real(k0), np.real(k1), np.real(k2))
    k_maj = (kl_max + k_vac) / 2

    def integrand_real(k_rho):
        kz = branch(1, k0, k_rho)
        k1z = branch(epsilon1, k0, k_rho)
        k2z = branch(epsilon2, k0, k_rho)

        rs_plus = rs(kz, k1z)
        rs_minus = rs(kz, k2z)
        rp_plus = rp(kz, k1z, epsilon1)
        rp_minus = rp(kz, k2z, epsilon2)

        # Integrand
        to_integrate = (
            1j
            / (4 * np.pi)
            * (
                Gs(kz, h, k_rho, rho, phi, rs_plus, rs_minus, rp_plus, rp_minus, z_ges, z_ab)
                - (kz**2 / k0**2) * Gp(kz, h, k_rho, rho, phi, rs_plus, rs_minus, rp_plus, rp_minus, z_ges, z_ab)
            )
            * (k_rho / kz)
            * np.exp(1j * kz * h)
        )

        return to_integrate

    # Real and imaginary part
    real_part_real = lambda k_rho: np.real(integrand_real(k_rho))
    imag_part_real = lambda k_rho: np.imag(integrand_real(k_rho))

    Re_real, _ = quad(real_part_real, 2 * k_maj, upper_limit * np.real(k0), limit=1000, epsrel=1e-9)
    Im_real, _ = quad(imag_part_real, 2 * k_maj, upper_limit * np.real(k0), limit=1000, epsrel=1e-9)

    return Re_real + 1j * Im_real


def green_tensor_scattered(
    r_a: np.ndarray, r_b: np.ndarray, omega: float, epsilon0: complex, epsilon1: complex, epsilon2: complex, h: float
) -> np.ndarray:
    """ Assemble the total scattering Green tensor

    Args:
        r_a: Position vector of atom A (m)
        r_b: Position vector of atom B (m)
        omega: Angular frequency in (1/s)
        epsilon1: Electric permittivity of the upper medium (dimensionless, complex)
        epsilon2: Electric permittivity of the lower medium (dimensionless, complex)
        h: Distance between the two surfaces (m)

    Returns: The 3x3 Scattering Green Tensor (general complex values) (1/m)
    """
    r = r_a - r_b
    rho = np.sqrt(r[0] ** 2 + r[1] ** 2)
    z_alpha = r_a[2]
    z_beta = r_b[2]
    z_ges = z_alpha + z_beta
    z_AB = r[2]
    phi = (
        np.arccos(r[0] / rho) if rho != 0 else 0
    )

    # Estimate the upper limit for the real axis integral
    k_vac = omega / const.c  # magnitude of wave vector in vacuum
    k0 = k_vac * np.sqrt(epsilon0 + 0j)
    upper_limit_real_integral = np.sqrt(
        (745 / (np.real(k0) * h)) ** 2 + 1
    )

    # Matrices of scattered Green's functions for s- and p-polarized light
    gs_matrix = np.array([[Gs_xx, Gs_xy, Gs_xz], [Gs_yx, Gs_yy, Gs_yz], [Gs_zx, Gs_zy, Gs_zz]])
    gp_matrix = np.array([[Gp_xx, Gp_xy, Gp_xz], [Gp_yx, Gp_yy, Gp_yz], [Gp_zx, Gp_zy, Gp_zz]])

    gt_total = np.zeros((3, 3), dtype=complex)

    for i in range(3):
        for j in range(3):
            gs_ij = gs_matrix[i][j]
            gp_ij = gp_matrix[i][j]
            g_ij_elliptic = elliptic_integral(
                omega, h, rho, phi, epsilon0, epsilon1, epsilon2, z_ges, z_AB, gs_ij, gp_ij
            )
            g_ij_real = real_axis_integral(
                omega, h, rho, phi, epsilon0, epsilon1, epsilon2, z_ges, z_AB, gs_ij, gp_ij, upper_limit_real_integral
            )
            gt_total[i][j] = g_ij_elliptic + g_ij_real

    return gt_total

def green_tensor_total(
    r_a: np.ndarray, r_b: np.ndarray, omega: float, epsilon0: complex, epsilon1: complex, epsilon2: complex, h: float
) -> np.ndarray:
    """ Assemble the total Green tensor
    Args:
        r_a: Position vector of atom A (m)
        r_b: Position vector of atom B (m)
        omega: Angular frequency in (1/s)
        epsilon1: Electric permittivity of the upper medium (dimensionless, complex)
        epsilon2: Electric permittivity of the lower medium (dimensionless, complex)
        h: Distance between the two surfaces (m)

    Returns: The 3x3 Total Green Tensor (general complex values) (1/m)
    """
    gt_scat = green_tensor_scattered(r_a, r_b, omega, epsilon0, epsilon1, epsilon2, h)
    gt_hom = green_tensor_homogeneous(r_a, r_b, omega, epsilon0)

    return gt_scat + gt_hom
