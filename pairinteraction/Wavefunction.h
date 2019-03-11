/*
 * Copyright (c) 2017 Sebastian Weber, Henri Menke. All rights reserved.
 *
 * This file is part of the pairinteraction library.
 *
 * The pairinteraction library is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * The pairinteraction library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with the pairinteraction library. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef WAVEFUNCTION_H
#define WAVEFUNCTION_H

#include "QuantumDefect.h"
#include "dtypes.h"
#include <string>
#include <vector>

// --- Numerov's method ---

namespace model_potential {
/** \brief Calculate the model potentials.
 *
 * This function calculates the model potentials and returns \f$
 * V_{\mathrm{mod}}(r) = V_{\mathrm{C}}(r) + V_{\mathrm{P}}(r) +
 * V_{\mathrm{s.o.}}(r) \f$
 *
 * \f[
 * \begin{aligned}
 *      V_{\mathrm{C}}(r) &= - \frac{e^2}{4\pi\varepsilon_0} \frac{1 +
 * (Z-1)\mathrm{e}^{-\alpha_1 r} - r(\alpha_3+\alpha_4 r)\mathrm{e}^{-\alpha_2
 * r}}{r} \\
 *      V_{\mathrm{P}}(r) &= -\frac{e^2}{(4\pi\varepsilon_0)^2}
 * \frac{\alpha_d}{2 r^4} \left[ 1 - \mathrm{e}^{-(r/r_c)^6} \right] \\
 *      V_{\mathrm{s.o.}}(r) &= \frac{1}{2} \left( \frac{e^2}{4\pi\varepsilon_0}
 * \right) \left( \frac{g_s}{2 m_e^2 c^2} \right)
 * \frac{\boldsymbol{l}\cdot\boldsymbol{s}}{r^3} \end{aligned} \f]
 *
 * \param[in] qd    Quantum defect data (parameters for the potentials)
 * \param[in] x     Radial position
 * \returns Model potentials evaluted at position \p x with parameters \p qd
 */
double V(QuantumDefect const &qd, double x);

/** \brief Helper function for %Numerov's method.s
 *
 * %Numerov's method solves \f$ y''(x) + g(x) y(x) = 0\f$. This function
 * implements \f$ g(x) \f$.
 *
 * \f[
 *      g(x) = \frac{(2 l + 1/2)^2}{x} + 8 x (V_{\mathrm{mod}}(x) - E)
 * \f]
 *
 * \param[in] qd    Quantum defect data (parameters for the potentials)
 * \param[in] x     Radial position
 * \returns Scaling coefficient evaluted at position \p x with parameters \p qd
 */
double g(QuantumDefect const &qd, double x);
} // namespace model_potential

/** \brief %Numerov's method
 *
 * This class implements %Numerov's method using a couple of helper functions.
 * %Numerov's method solves a differential equation of the form
 * \f[
 *      y''(x) + g(x)y(x) = 0
 * \f]
 * The equation is solved by iteratively computing
 * \f[
 *      y_{n+1} = \frac{\left( 2 + \frac{5 h^2}{6} g_n \right) y_n - \left( 1 +
 * \frac{h^2}{12} g_{n-1} \right) y_{n-1}}{\left( 1 + \frac{h^2}{12} g_n
 * \right)} \f]
 *
 * with \f$ y_n = y(x_n) \f$ and \f$ g_n = g(x_n) \f$. Clearly this equation
 * needs two initial conditions \f$ y_0 \f$ and \f$ y_1 \f$. We integrate the
 * equation from outside to inside. Since the wavefunction hs to decay to zero
 * at infinity we choose \f$ y_0 = 0 \f$. Now depending on the number of knots
 * of wavefunction we decide whether to set \f$ y_1 = \pm\varepsilon \f$.
 */
class Numerov {
    QuantumDefect const &qd;
    eigen_dense_double_t xy;

public:
    /** \brief Integration step size */
    static constexpr double const dx = 0.01;

    /** \brief Constructor
     *
     * The constructor initializes the \p qd member and sets the initial
     * condition for the integrator. It determines the outer and inner
     * integration bounds by semiclassical and empirical arguments.
     *
     * \param[in] qd    Quantum defect data (parameters for the potentials)
     */
    Numerov(QuantumDefect const &qd);

    /** \brief Perform the integration
     *
     * This performs the integration using %Numerov's method.
     *
     * \returns vector with wavefunction amplitude
     */
    eigen_dense_double_t integrate();

    /** \brief Power kernel for matrix elements
     *
     * The power kernel accounts for the fact that in %Numerov's method the
     * domain is square root scaled. This is important for the calculation of
     * matrix elements.
     *
     * \param[in] power     exponent of r
     * \returns power kernel
     */
    constexpr static inline int power_kernel(int power) { return 2 * power + 2; }
};

// --- Whittaker method ---

namespace whittaker_functions {
/** \brief Compute the confluent hypergeometric function
 *
 * This is merely a wrapper around
 * <a
 * href="https://www.gnu.org/software/gsl/manual/html_node/Hypergeometric-Functions.html">
 * <tt>gsl_sf_hyperg_U</tt>
 * </a>
 *
 * \param[in] a     see documentation for <tt>gsl_sf_hyperg_U</tt>
 * \param[in] b     see documentation for <tt>gsl_sf_hyperg_U</tt>
 * \param[in] z     see documentation for <tt>gsl_sf_hyperg_U</tt>
 * \returns U(a,b,z)
 */
double HypergeometricU(double a, double b, double z);

/** \brief Compute the %Whittaker function
 *
 * The %Whittaker function is defined in terms of the confluent hypergeometric
 * function.
 * \f[
 *        W_{k,m}(z) = \mathrm{e}^{-z/2} z^{m+1/2} U\left(m-k+\frac{1}{2}, 1+2m,
 * z\right) \f]
 *
 * \param[in] k     parameter k from %Whittaker's equation
 * \param[in] m     parameter m from %Whittaker's equation
 * \param[in] z     radial position
 * \returns W(k,m,z)
 */
double WhittakerW(double k, double m, double z);

/** \brief Radial wavefunctions from %Whittaker's function
 *
 * The radial wavefunction of the Schrödinger equation with the Coulomb
 * potential can be expressed in terms of the %Whittaker function.
 * \f[
 *        R_{\nu,l}(r) = \frac{1}{\sqrt{\nu^2 \Gamma(\nu+l+1) \Gamma(\nu-l)}}
 * W_{\nu, l+1/2}\left(\frac{2r}{\nu}\right) \f]
 *
 * \param[in] nu    effective principal quantum number \f$ \nu = n -
 * \delta_{n,l} \f$ \param[in] l     angular quantum number \param[in] r radial
 * position \returns R(nu,l,r)
 */
double RadialWFWhittaker(double r, double nu, int l);
} // namespace whittaker_functions

class Whittaker {
    QuantumDefect const &qd;
    eigen_dense_double_t xy;

public:
    /** \brief Integration step size */
    static constexpr double const dx = 0.01;

    /** \brief Constructor
     *
     * The constructor initializes the \p qd member and sets the initial
     * condition for the integrator. It determines the outer integration bound
     * by semiclassical and empirical arguments. The inner intergration
     * boundary is set to 1.
     *
     * \param[in] qd    Quantum defect data (parameters for the potentials)
     */
    Whittaker(QuantumDefect const &qd);

    /** \brief Evaluate the radial wavefunction
     *
     * This does not perform an integration but merely evaluates the %Whittaker
     * functions over the domain.
     *
     * \returns vector with wavefunction amplitude
     */
    eigen_dense_double_t integrate();

    /** \brief Power kernel for matrix elements
     *
     * The power kernel accounts for the fact that for the %Whittaker functions
     * the domain is linear scaled. This is important for the calculation of
     * matrix elements.
     *
     * \param[in] power     exponent of r
     * \returns power kernel
     */
    constexpr static inline double power_kernel(int power) { return 2 * power + 1; }
};

// --- Matrix element calculation ---

/** \brief Find and return index
 *
 * Find a value in an Eigen matrix and return its index in the
 * container to the caller. This is intended to work with
 * Eigen::Vector or columns of Eigen::Matrix.
 *
 * \param[in] x     An Eigen matrix type
 * \param[in] d     value to find
 * \returns index of d
 * \throws std::runtime_error if the value can't be found
 */
template <typename T>
int findidx(T const &x, typename T::Scalar const &d) {
    int L = 0;
    int R = x.rows() - 1;
    for (;;) {
        if (L > R) {
            throw std::runtime_error("Search failed");
        }
        int m = (L + R) / 2;
        if (x(m) < d) {
            L = m + 1;
        }
        if (x(m) > d) {
            R = m - 1;
        }
        if (x(m) == d) {
            return m;
        }
    }
}

/** \brief Compute radial matrix elements
 *
 * The radial matrix element can be calculated from the integral
 * \f[
 *      \langle nlj| \hat{p}^{\mathrm{rad}}_\kappa |n'l'j' \rangle =
 *      \int \Psi^{\text{rad}}_{nlj}(r) \Psi^{\text{rad}}_{n'l'j'}(r)
 * r^{2+\kappa} \; dr \f] The part \f$ r^{2+\kappa} \f$ varies with rescaling
 * the domain and is this given by the power_kernel function.
 *
 * \tparam T        Method for calculating the wavefunction
 * \param[in] qd1   Quantum  defect data for first atom
 * \param[in] power Exponent kappa
 * \param[in] qd2   Quantum  defect data for second atom
 * \returns Radial matrix element
 */
template <typename T>
double IntegrateRadialElement(QuantumDefect const &qd1, int power, QuantumDefect const &qd2) {
    T N1(qd1);
    T N2(qd2);

    auto const &xy1 = N1.integrate();
    auto const &xy2 = N2.integrate();
    auto const dx = T::dx;

    auto const xmin = xy1(0, 0) >= xy2(0, 0) ? xy1(0, 0) : xy2(0, 0);
    auto const xmax = xy1(xy1.rows() - 1, 0) <= xy2(xy2.rows() - 1, 0) ? xy1(xy1.rows() - 1, 0)
                                                                       : xy2(xy2.rows() - 1, 0);

    double mu = 0;
    // If there is an overlap, calculate the matrix element
    if (xmin <= xmax) {
        int start1 = findidx(xy1.col(0), xmin);
        int end1 = findidx(xy1.col(0), xmax);
        int start2 = findidx(xy2.col(0), xmin);
        int end2 = findidx(xy2.col(0), xmax);

        int i1, i2;
        for (i1 = start1, i2 = start2; i1 < end1 && i2 < end2; ++i1, ++i2) {
            mu += xy1(i1, 1) * xy2(i2, 1) * std::pow(xy1(i1, 0), T::power_kernel(power)) * dx;
        }
        mu = 2 * mu;
    }

    // The radial matrix element is returned in atomic units
    return mu;
}

#endif // WAVEFUNCTION_H
