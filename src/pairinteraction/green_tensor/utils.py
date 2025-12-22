# SPDX-FileCopyrightText: 2024 PairInteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

import numpy as np
import scipy.constants as const
from numpy.typing import NDArray
from scipy.special import jn
from scipy.integrate import quad


# TODO implement utils for green tensors

# ---- Homohgeneous Green Tensor ----
''' The function used is from equation 2 of the paper: "Dispersionless subradiant photon storage in one-dimensional emitter chains" 

    Args: 
        r_A: Position vector of atom A in meters
        r_B: Position vector of atom B in meters
        omega: Angular frequency in 1/seconds
        
    Returns: The 3x3 Homogeneous Green Tensor in cartesian coordinates with in general complex values 
    '''
# Homogeneous Green Tensor for two atoms in free space in cartesian coordinates
def green_tensor_homogeneous(r_A: NDArray, r_B: NDArray, omega: float, epsilon0: complex) -> NDArray: 
    k_vac = omega / const.c # magnitude of wave vector in vacuum
    k0 = k_vac * np.sqrt(epsilon0 + 0j) # magnitude of wave vector in medium with permittivity epsilon0
    r = r_A - r_B # vector distance between the two atoms
    r_norm = np.linalg.norm(r)    

    prefactor = np.exp(1j * k0 * r_norm) / (4 * np.pi * k0**2 * r_norm**3)
    term1 = (k0**2 * r_norm**2 + 1j * k0 * r_norm - 1) * np.eye(3)
    term2 = (-k0**2 * r_norm**2 - 3j * k0 * r_norm + 3) * np.outer(r, r) / r_norm**2
    G_hom = prefactor * (term1 + term2)

    return G_hom

# ---- Scattering Green Tensor ----
''' The functions used are from Appendix B of the paper: "Modified dipole-dipole interaction and dissipation in an atomic ensemble near surfaces" '''

# ---- Functions needed for the calculation of the matrix elements of the scattering Green Tensor ----
''' In general the scattering Green Tensor is a 3x3 matrix with each element being either zero or an integral over the in-plane wave vector k_rho. 
    The scattering Green Tensor can be devided into contributions from s- and p-polarized light. These contributions themselves are 3x3 matrices.

    Some of the matrix elements share common prefactors and terms, which are calculated in the following functions.'''

# Branch of the wave vector component perpendicular to the surfaces
''' This function ensures that the imaginary part of the perpendicular wave vector component is always positive.
    From my experience the function np.sqrt already does this correctly, but to be sure I implemented this function. 

    Args: 
        epsilon: Electric permittivity of the medium (dimensionless, complex)
        k: Magnitude of wave number in vacuum (1/m)
        k_rho: In-plane wave vector component (1/m)
    '''
def branch(epsilon, k, k_rho):
    val = np.sqrt(epsilon * k**2 - k_rho**2 + 0j)
    return val
#    return val if np.imag(val) >= 0 else -val
#    return abs(np.real(val)) + 1j * abs(np.imag(val))

# Fresnel reflection coefficient for s-polarized light
def rs(kz: complex, k1z: complex) -> complex:
    return (kz - k1z) / (kz + k1z)

# Fresnel reflection coefficient for p-polarized light
def rp(kz: complex, k1z: complex, epsilon: complex) -> complex:
    return (epsilon * kz - k1z) / (epsilon * kz + k1z)

# Prefactors for matrix elements of scattering Green's Tensor
def D(r_plus: complex, r_minus: complex, kz: complex, h: float) -> complex:
    return 1 - r_plus * r_minus * np.exp(2j * kz * h)

def A_plus(r_plus: complex, r_minus: complex, kz: complex, h: float, z_ges: float, z_AB: float) -> complex:
    return (r_minus * np.exp(1j * kz * (z_ges - h)) + r_plus * np.exp(-1j * kz * (z_ges - h)) + 2 * r_plus * r_minus * np.cos(kz * z_AB) * np.exp(1j * kz * h)) / D(r_plus, r_minus, kz, h)

def A_minus(r_plus: complex, r_minus: complex, kz: complex, h: float, z_ges: float, z_AB: float) -> complex:
    return (r_minus * np.exp(1j * kz * (z_ges - h)) + r_plus * np.exp(-1j * kz * (z_ges - h)) - 2 * r_plus * r_minus * np.cos(kz * z_AB) * np.exp(1j * kz * h)) / D(r_plus, r_minus, kz, h)

def B_plus(r_plus: complex, r_minus: complex, kz: complex, h: float, z_ges: float, z_AB: float) -> complex:
    return (r_minus * np.exp(1j * kz * (z_ges - h)) + r_plus * np.exp(-1j * kz * (z_ges - h)) + 2j * r_plus * r_minus * np.sin(kz * z_AB) * np.exp(1j * kz * h)) / D(r_plus, r_minus, kz, h)

def B_minus(r_plus: complex, r_minus: complex, kz: complex, h: float, z_ges: float, z_AB: float) -> complex:
    return (r_minus * np.exp(1j * kz * (z_ges - h)) + r_plus * np.exp(-1j * kz * (z_ges - h)) - 2j * r_plus * r_minus * np.sin(kz * z_AB) * np.exp(1j * kz * h)) / D(r_plus, r_minus, kz, h)

# ---- Functions to calculate individual matrix elements of the scattering Green Tensor ----

# Diagonal matrix elements of scattering Green Tensor
def Gs_xx(kz: complex, h: float, k_rho: complex, rho: float, phi: float, rs_plus: complex, rs_minus: complex, rp_plus: complex, rp_minus: complex, z_ges: float, z_AB: float) -> complex:
    # Bessel functions
    J0 = jn(0, k_rho * rho)
    J2 = jn(2, k_rho * rho)

    As_plus = A_plus(rs_plus, rs_minus, kz, h, z_ges, z_AB)

    return As_plus / 2 * (J0 + J2 * np.cos(2 * phi))

def Gs_yy(kz: complex, h: float, k_rho: complex, rho: float, phi: float, rs_plus: complex, rs_minus: complex, rp_plus: complex, rp_minus: complex, z_ges: float, z_AB: float) -> complex:
    # Bessel functions
    J0 = jn(0, k_rho * rho)
    J2 = jn(2, k_rho * rho)

    As_plus = A_plus(rs_plus, rs_minus, kz, h, z_ges, z_AB)

    return As_plus / 2 * (J0 - J2 * np.cos(2 * phi))

def Gs_zz(kz: complex, h: float, k_rho: complex, rho: float, phi: float, rs_plus: complex, rs_minus: complex, rp_plus: complex, rp_minus: complex, z_ges: float, z_AB: float) -> complex:
    return 0

def Gp_xx(kz: complex, h: float, k_rho: complex, rho: float, phi: float, rs_plus: complex, rs_minus: complex, rp_plus: complex, rp_minus: complex, z_ges: float, z_AB: float) -> complex:
    # Bessel functions
    J0 = jn(0, k_rho * rho)
    J2 = jn(2, k_rho * rho)

    Ap_minus = A_minus(rp_plus, rp_minus, kz, h, z_ges, z_AB)

    return Ap_minus / 2 * (J0 - J2 * np.cos(2 * phi))

def Gp_yy(kz: complex, h: float, k_rho: complex, rho: float, phi: float, rs_plus: complex, rs_minus: complex, rp_plus: complex, rp_minus: complex, z_ges: float, z_AB: float) -> complex:
    # Bessel functions
    J0 = jn(0, k_rho * rho)
    J2 = jn(2, k_rho * rho)

    Ap_minus = A_minus(rp_plus, rp_minus, kz, h, z_ges, z_AB)

    return Ap_minus / 2 * (J0 + J2 * np.cos(2 * phi))

def Gp_zz(kz: complex, h: float, k_rho: complex, rho: float, phi: float, rs_plus: complex, rs_minus: complex, rp_plus: complex, rp_minus: complex, z_ges: float, z_AB: float) -> complex:
    # Bessel function
    J0 = jn(0, k_rho * rho)

    Ap_plus = A_plus(rp_plus, rp_minus, kz, h, z_ges, z_AB)

    return - (k_rho**2 / kz**2) * Ap_plus * J0

# Off-diagonal elements of scattered Green's Tensor
def Gs_xy(kz: complex, h: float, k_rho: complex, rho: float, phi: float, rs_plus: complex, rs_minus: complex, rp_plus: complex, rp_minus: complex, z_ges: float, z_AB: float) -> complex:
    J2 = jn(2, k_rho * rho)
    As_plus = A_plus(rs_plus, rs_minus, kz, h, z_ges, z_AB)

    return - As_plus / 2 * J2 * np.sin(2 * phi)

def Gs_yx(kz: complex, h: float, k_rho: complex, rho: float, phi: float, rs_plus: complex, rs_minus: complex, rp_plus: complex, rp_minus: complex, z_ges: float, z_AB: float) -> complex:
    J2 = jn(2, k_rho * rho)
    As_plus = A_plus(rs_plus, rs_minus, kz, h, z_ges, z_AB)

    return - As_plus / 2 * J2 * np.sin(2 * phi)

def Gs_xz(kz: complex, h: float, k_rho: complex, rho: float, phi: float, rs_plus: complex, rs_minus: complex, rp_plus: complex, rp_minus: complex, z_ges: float, z_AB: float) -> complex:
    return 0

def Gs_zx(kz: complex, h: float, k_rho: complex, rho: float, phi: float, rs_plus: complex, rs_minus: complex, rp_plus: complex, rp_minus: complex, z_ges: float, z_AB: float) -> complex:
    return 0

def Gs_yz(kz: complex, h: float, k_rho: complex, rho: float, phi: float, rs_plus: complex, rs_minus: complex, rp_plus: complex, rp_minus: complex, z_ges: float, z_AB: float) -> complex:
    return 0

def Gs_zy(kz: complex, h: float, k_rho: complex, rho: float, phi: float, rs_plus: complex, rs_minus: complex, rp_plus: complex, rp_minus: complex, z_ges: float, z_AB: float) -> complex:
    return 0

def Gp_xy(kz: complex, h: float, k_rho: complex, rho: float, phi: float, rs_plus: complex, rs_minus: complex, rp_plus: complex, rp_minus: complex, z_ges: float, z_AB: float) -> complex:
    J2 = jn(2, k_rho * rho)
    Ap_minus = A_minus(rp_plus, rp_minus, kz, h, z_ges, z_AB)

    return Ap_minus / 2 * J2 * np.sin(2 * phi)

def Gp_xz(kz: complex, h: float, k_rho: complex, rho: float, phi: float, rs_plus: complex, rs_minus: complex, rp_plus: complex, rp_minus: complex, z_ges: float, z_AB: float) -> complex:
    J1 = jn(1, k_rho * rho)
    Bp_plus = B_plus(rp_plus, rp_minus, kz, h, z_ges, z_AB)

    return 1j * (k_rho / kz) * Bp_plus * J1 * np.cos(phi)

def Gp_yx(kz: complex, h: float, k_rho: complex, rho: float, phi: float, rs_plus: complex, rs_minus: complex, rp_plus: complex, rp_minus: complex, z_ges: float, z_AB: float) -> complex:
    J2 = jn(2, k_rho * rho)
    Ap_minus = A_minus(rp_plus,rp_minus, kz, h, z_ges, z_AB)

    return Ap_minus / 2 * J2 * np.sin(2 * phi)

def Gp_yz(kz: complex, h: float, k_rho: complex, rho: float, phi: float, rs_plus: complex, rs_minus: complex, rp_plus: complex, rp_minus: complex, z_ges: float, z_AB: float) -> complex:
    J1 = jn(1, k_rho * rho)
    Bp_plus = B_plus(rp_plus, rp_minus, kz, h, z_ges, z_AB)

    return 1j * (k_rho / kz) * Bp_plus * J1 * np.sin(phi)

def Gp_zx(kz: complex, h: float, k_rho: complex, rho: float, phi: float, rs_plus: complex, rs_minus: complex, rp_plus: complex, rp_minus: complex, z_ges: float, z_AB: float) -> complex:
    J1 = jn(1, k_rho * rho)
    Bp_minus = B_minus(rp_plus, rp_minus, kz, h, z_ges, z_AB)

    return -1j * (k_rho / kz) * Bp_minus * J1 * np.cos(phi)

def Gp_zy(kz: complex, h: float, k_rho: complex, rho: float, phi: float, rs_plus: complex, rs_minus: complex, rp_plus: complex, rp_minus: complex, z_ges: float, z_AB: float) -> complex:
    J1 = jn(1, k_rho * rho)
    Bp_plus = B_plus(rp_plus, rp_minus, kz, h, z_ges, z_AB)

    return -1j * (k_rho / kz) * Bp_plus * J1 * np.sin(phi)

# Matrices of scattered Green's functions for s- and p-polarized light
Gs_matrix = np.array([
            [Gs_xx, Gs_xy, Gs_xz], 
            [Gs_yx, Gs_yy, Gs_yz], 
            [Gs_zx, Gs_zy, Gs_zz]])

Gp_matrix = np.array([
            [Gp_xx, Gp_xy, Gp_xz], 
            [Gp_yx, Gp_yy, Gp_yz], 
            [Gp_zx, Gp_zy, Gp_zz]])

# ---- Functions to evaluate the integrals for the scattering Green Tensor ----
''' The integrals for the scattering Green Tensor are evaluated in two parts: 
    The first part is along an elliptical path from 0 to 2*k_maj in the complex plane to avoid singularities,
    and the second part is along the real axis from 2*k_maj to an upper limit. 
    
    The methods for evaluating these integrals are explained in the papers:
    "Accurate and efficient computation of the Green’s tensor for stratified media" section III.A
    "Challenges in Computational Electromagnetics: Analysis and Optimization of Planar Multilayered Structures" section 2.4.1
    '''

# Function to evaluate the elliptic part of the integral
''' Args:
        omega: Angular frequency in 1/seconds
        h: Distance between the two surfaces in meters
        rho: In-plane distance between the two atoms in meters
        phi: Angle between the in-plane distance vector and the x-axis in radians
        epsilon0: Electric permittivity of the medium between the two surfaces (dimensionless, complex)
        epsilon1: Electric permittivity of the upper medium (dimensionless, complex)
        epsilon2: Electric permittivity of the lower medium (dimensionless, complex)
        z_ges: Sum of the z-positions of the two atoms in meters
        z_AB: Difference of the z-positions of the two atoms in meters
        Gs: Function to calculate the s-polarized Green Tensor matrix elements
        Gp: Function to calculate the p-polarized Green Tensor matrix elements

    Returns: The value of the integral along the elliptical path as a complex number
    '''
def elliptic_integral(omega: float, h: float, rho: float, phi: float, epsilon0: complex, epsilon1: complex, epsilon2: complex, z_ges: float, z_AB: float, Gs, Gp) -> complex:
    k_vac = omega / const.c # magnitude of wave vector in vacuum
    k0 = k_vac * np.sqrt(epsilon0 + 0j)

    # Elliptical path in complex plane to avoid singularities (Integral from 0 to 2k_maj)
    k1 = k_vac * np.sqrt(epsilon1 + 0j)
    k2 = k_vac * np.sqrt(epsilon2 + 0j)
    kl_max = max(np.real(k0), np.real(k1), np.real(k2))

    k_maj = (kl_max + k_vac) / 2 # major axis of ellipse
    k_min = min(k_vac, 1 / rho) if rho != 0 else k_vac # minor axis of ellipse, if clause only needed if systems with only one atom is considered as well

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
        to_integrate = 1j / (4 * np.pi) * (Gs(kz, h, k_rho, rho, phi, rs_plus, rs_minus, rp_plus, rp_minus, z_ges, z_AB) - (kz**2 / k0**2) * Gp(kz, h, k_rho, rho, phi, rs_plus, rs_minus, rp_plus, rp_minus, z_ges, z_AB)) * (k_rho / kz) * np.exp(1j * kz * h) * dk_rho

        return to_integrate
    
    # Real and imaginary part
    real_part_ellipse = lambda t: np.real(integrand_ellipse(t))
    imag_part_ellipse = lambda t: np.imag(integrand_ellipse(t))

    Re_ellipse, _ = quad(real_part_ellipse, np.pi, 0, epsrel=1e-9, limit=1000)
    Im_ellipse, _ = quad(imag_part_ellipse, np.pi, 0, epsrel=1e-9, limit=1000)

    return (Re_ellipse + 1j * Im_ellipse)

# Function to evaluate the integral along the real axis
''' Args:
        same as for elliptic_integral plus:
        upper_limit: Upper limit of the integral along the real axis (in units of k)'''
def real_axis_integral(omega: float, h: float, rho: float, phi: float, epsilon0: complex, epsilon1: complex, epsilon2: complex, z_ges: float, z_AB: float, Gs, Gp, upper_limit: float) -> complex:
    k_vac = omega / const.c # magnitude of wave vector in vacuum
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
        to_integrate = 1j / (4 * np.pi) * (Gs(kz, h, k_rho, rho, phi, rs_plus, rs_minus, rp_plus, rp_minus, z_ges, z_AB) - (kz**2 / k0**2) * Gp(kz, h, k_rho, rho, phi, rs_plus, rs_minus, rp_plus, rp_minus, z_ges, z_AB)) * (k_rho / kz) * np.exp(1j * kz * h)

        return to_integrate
    
    # Real and imaginary part
    real_part_real = lambda k_rho: np.real(integrand_real(k_rho))
    imag_part_real = lambda k_rho: np.imag(integrand_real(k_rho))

    Re_real, _ = quad(real_part_real, 2 * k_maj, upper_limit * np.real(k0), limit=1000, epsrel=1e-9)
    Im_real, _ = quad(imag_part_real, 2 * k_maj, upper_limit * np.real(k0), limit=1000, epsrel=1e-9)

    return (Re_real + 1j * Im_real)

# Assemble the total scattering Green Tensor
''' Args:
        r_A: Position vector of atom A in meters
        r_B: Position vector of atom B in meters
        omega: Angular frequency in 1/seconds
        epsilon1: Electric permittivity of the upper medium (dimensionless, complex)
        epsilon2: Electric permittivity of the lower medium (dimensionless, complex)
        h: Distance between the two surfaces in meters

    Returns: The 3x3 Scattering Green Tensor in cartesian coordinates with in general complex values 
    '''
def green_tensor_scattered(r_A: np.ndarray, r_B: np.ndarray, omega: float, epsilon0: complex, epsilon1: complex, epsilon2: complex, h: float) -> np.ndarray:
    r = r_A - r_B # vector distance between the two atoms
    rho = np.sqrt(r[0]**2 + r[1]**2) # in-plane distance between the two atoms
    z_alpha = r_A[2] # z-coordinate of atom A
    z_beta = r_B[2] # z-coordinate of atom B
    z_ges = z_alpha + z_beta # sum of the z-positions of the two atoms
    z_AB = r[2] # difference of the z-positions of the two atoms
    phi = np.arccos(r[0] / rho) if rho != 0 else 0 # angle between the in-plane distance vector and the x-axis, if clause only needed if systems with only one atom is considered as well
    
    # Estimate the upper limit for the real axis integral
    k_vac = omega / const.c # magnitude of wave vector in vacuum
    k0 = k_vac * np.sqrt(epsilon0 + 0j)
    upper_limit_real_integral = np.sqrt((745 / (np.real(k0) * h))**2 + 1)   # upper limit for the real axis integral, evaluated such that the leading term of the integrand (the exponential function,
                                                                            # exp(1j * kz * h) = exp(-sqrt(k_rho**2 - k**2) * h) for large k_rho) becomes zero at the upper limit
                                                                            # For simplicity, the upper limit is evaluated using only the real part of k0.

    G_total = np.zeros((3, 3), dtype=complex)

    for i in range(3):
        for j in range(3):
            Gs_ij = Gs_matrix[i][j]
            Gp_ij = Gp_matrix[i][j]
            G_ij_elliptic = elliptic_integral(omega, h, rho, phi, epsilon0, epsilon1, epsilon2, z_ges, z_AB, Gs_ij, Gp_ij)
            G_ij_real = real_axis_integral(omega, h, rho, phi, epsilon0, epsilon1, epsilon2, z_ges, z_AB, Gs_ij, Gp_ij, upper_limit_real_integral)
            G_total[i][j] = G_ij_elliptic + G_ij_real

    return G_total

# ---- Assemble Green Tensor ----
''' The total Green Tensor is the sum of the homogeneous and scattering Green Tensors '''
# Total Green Tensor for two atoms near one or two surfaces in cartesian coordinates
def green_tensor_total(r_A: np.ndarray, r_B: np.ndarray, omega: float, epsilon0: complex, epsilon1: complex, epsilon2: complex, h: float) -> np.ndarray:
    G_scat = green_tensor_scattered(r_A, r_B, omega, epsilon0, epsilon1, epsilon2, h)
    G_hom = green_tensor_homogeneous(r_A, r_B, omega, epsilon0)

    return G_scat + G_hom

''' For completeness, these two functions are included here, although they are not used in the current implementation of the Green Tensor classes,
    since they are only needed, if only one atom in the system is considered. '''

# Homogeneous Green Tensor for one atom in cartesian coordinates
def green_tensor_homogeneous_one_atom(omega, epsilon0): 
    k_vac = omega / const.c
    k0 = k_vac * np.sqrt(epsilon0 + 0j)

    return 1j * k0 / (6 * np.pi) * np.eye(3)

# Total Green Tensor for one atom near one or two surfaces in cartesian coordinates
def green_tensor_total_one_atom(r_A, omega, epsilon0, epsilon1, epsilon2, h):
    G_scat = green_tensor_scattered(r_A, r_A, omega, epsilon0, epsilon1, epsilon2, h)
    G_hom = green_tensor_homogeneous_one_atom(omega, epsilon0)

    return G_scat + G_hom

''' Findings: 
            It seems, that if the height between the two surfaces is large (e.g. several wavelengths), the oscillarity behaviour of the integrand
            becomes less problematic, since they mainly occur on the elliptic path and the integral evaluation is faster.
            If the height is small (e.g. much smaller than a wavelength), the oscillations extend far into the real axis integral and the evaluation
            becomes slow. Possible improvements could be to evaluate the integral along the real axis with a different method, e.g. using
            Extrapolation methods or DE rules, which I failed to make work yet.
            '''