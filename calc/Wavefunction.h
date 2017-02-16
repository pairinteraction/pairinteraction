/*
 * Copyright (c) 2017 Sebastian Weber, Henri Menke. All rights reserved.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#ifndef WAVEFUNCTION_H
#define WAVEFUNCTION_H

#include <string>
#include <vector>
#include "dtypes.h"
#include "QuantumDefect.h"

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
     *      V_{\mathrm{C}}(r) &= - \frac{e^2}{4\pi\varepsilon_0} \frac{1 + (Z-1)\mathrm{e}^{-\alpha_1 r} - r(\alpha_3+\alpha_4 r)\mathrm{e}^{-\alpha_2 r}}{r} \\
     *      V_{\mathrm{P}}(r) &= -\frac{e^2}{(4\pi\varepsilon_0)^2} \frac{\alpha_d}{2 r^4} \left[ 1 - \mathrm{e}^{-(r/r_c)^6} \right] \\
     *      V_{\mathrm{s.o.}}(r) &= \frac{1}{2} \left( \frac{e^2}{4\pi\varepsilon_0} \right) \left( \frac{g_s}{2 m_e^2 c^2} \right) \frac{\boldsymbol{l}\cdot\boldsymbol{s}}{r^3}
     * \end{aligned}
     * \f]
     *
     * \param[in] qd    Quantum defect data (parameters for the potentials)
     * \param[in] x     Radial position
     * \returns Model potentials evaluted at position \p x with parameters \p qd
     */
    real_t V(QuantumDefect const& qd, real_t x);

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
    real_t g(QuantumDefect const& qd, real_t x);
}

/** \brief %Numerov's method
 *
 * This class implements %Numerov's method using a couple of helper functions.
 * %Numerov's method solves a differential equation of the form
 * \f[
 *      y''(x) + g(x)y(x) = 0
 * \f]
 * The equation is solved by iteratively computing
 * \f[
 *      y_{n+1} = \frac{\left( 2 + \frac{5 h^2}{6} g_n \right) y_n - \left( 1 + \frac{h^2}{12} g_{n-1} \right) y_{n-1}}{\left( 1 + \frac{h^2}{12} g_n \right)}
 * \f]
 *
 * with \f$ y_n = y(x_n) \f$ and \f$ g_n = g(x_n) \f$. Clearly this equation
 * needs two initial conditions \f$ y_0 \f$ and \f$ y_1 \f$. We integrate the
 * equation from outside to inside. Since the wavefunction hs to decay to zero
 * at infinity we choose \f$ y_0 = 0 \f$. Now depending on the number of knots
 * of wavefunction we decide whether to set \f$ y_1 = \pm\varepsilon \f$.
 */
class Numerov {
    QuantumDefect const& qd;
    std::vector<real_t> x;
public:
    /** \brief Integration step size */
    real_t const dx;

    /** \brief Constructor
     *
     * The constructor initializes the \p qd member and sets the initial
     * condition for the integrator. It determines the outer and inner
     * integration bounds by semiclassical and empirical arguments.
     *
     * \param[in] qd    Quantum defect data (parameters for the potentials)
     */
    Numerov(QuantumDefect const& qd);

    /** \brief x values
     *
     * Returns the domain on which the wavefunction was evaluated.
     *
     * \returns vector with radial points
     */
    std::vector<real_t> axis() const;

    /** \brief Perform the integration
     *
     * This performs the integration using %Numerov's method.
     *
     * \returns vector with wavefunction amplitude
     */
    std::vector<real_t> integrate();

    /** \brief Power kernel for matrix elements
     *
     * The power kernel accounts for the fact that in %Numerov's method the
     * domain is square root scaled. This is important for the calculation of
     * matrix elements.
     *
     * \param[in] power     exponent of r
     * \returns power kernel
     */
    constexpr static inline int power_kernel(int power)
    {
        return 2*power+2;
    }
};

// --- Whittaker method ---

namespace whittaker_functions {
    /** \brief Compute the confluent hypergeometric function
     *
     * This is merely a wrapper around
     * <a href="https://www.gnu.org/software/gsl/manual/html_node/Hypergeometric-Functions.html">
     * <tt>gsl_sf_hyperg_U</tt>
     * </a>
     *
     * \param[in] a     see documentation for <tt>gsl_sf_hyperg_U</tt>
     * \param[in] b     see documentation for <tt>gsl_sf_hyperg_U</tt>
     * \param[in] z     see documentation for <tt>gsl_sf_hyperg_U</tt>
     * \returns U(a,b,z)
     */
    real_t HypergeometricU(real_t a, real_t b, real_t z);

    /** \brief Compute the %Whittaker function
     *
     * The %Whittaker function is defined in terms of the confluent hypergeometric
     * function.
     * \f[
     *        W_{k,m}(z) = \mathrm{e}^{-z/2} z^{m+1/2} U\left(m-k+\frac{1}{2}, 1+2m, z\right)
     * \f]
     *
     * \param[in] k     parameter k from %Whittaker's equation
     * \param[in] m     parameter m from %Whittaker's equation
     * \param[in] z     radial position
     * \returns W(k,m,z)
     */
    real_t WhittakerW(real_t k, real_t m, real_t z);

    /** \brief Radial wavefunctions from %Whittaker's function
     *
     * The radial wavefunction of the Schrödinger equation with the Coulomb
     * potential can be expressed in terms of the %Whittaker function.
     * \f[
     *        R_{\nu,l}(r) = (\nu^2 \Gamma(\nu+l+1) \Gamma(\nu-l))^{-1/2} W_{\nu, l+1/2}\left(\frac{2r}{\nu}\right)
     * \f]
     *
     * \param[in] nu    effective principal quantum number \f$ \nu = n - \delta_{n,l} \f$
     * \param[in] l     angular quantum number
     * \param[in] r     radial position
     * \returns R(nu,l,r)
     */
    real_t RadialWFWhittaker(real_t r, real_t nu, int l);
}

class Whittaker {
    QuantumDefect const& qd;
    std::vector<real_t> x;
public:
    /** \brief Integration step size */
    real_t const dx;

    /** \brief Constructor
     *
     * The constructor initializes the \p qd member and sets the initial
     * condition for the integrator. It determines the outer integration bound
     * by semiclassical and empirical arguments. The inner intergration
     * boundary is set to 1.
     *
     * \param[in] qd    Quantum defect data (parameters for the potentials)
     */
    Whittaker(QuantumDefect const& qd);

    /** \brief x values
     *
     * Returns the domain on which the wavefunction was evaluated.
     *
     * \returns vector with radial points
     */
    std::vector<real_t> axis() const;

    /** \brief Evaluate the radial wavefunction
     *
     * This does not perform an integration but merely evaluates the %Whittaker
     * functions over the domain.
     *
     * \returns vector with wavefunction amplitude
     */
    std::vector<real_t> integrate();

    /** \brief Power kernel for matrix elements
     *
     * The power kernel accounts for the fact that for the %Whittaker functions
     * the domain is linear scaled. This is important for the calculation of
     * matrix elements.
     *
     * \param[in] power     exponent of r
     * \returns power kernel
     */
    constexpr static inline real_t power_kernel(int power)
    {
        return 1.5*power;
    }
};


// --- Matrix element calculation ---


/** \brief Find and return index
 *
 * Find a value in a random access container and return its index in the
 * container to the caller. This is intended to work with std::vector.
 *
 * \param[in] x     random access container
 * \param[in] d     value to find
 * \returns index of d or x.size()
 */
template < typename T >
size_t findidx(T const& x, typename T::value_type const& d) {
    auto it = std::find(std::begin(x), std::end(x), d);
    return std::distance(std::begin(x), it);
}


/** \brief Compute radial matrix elements
 *
 * The radial matrix element can be calculated from the integral
 * \f[
 *      \langle nlj| \hat{p}^{\mathrm{rad}}_\kappa |n'l'j' \rangle =
 *      \int \Psi^{\text{rad}}_{nlj}(r) \Psi^{\text{rad}}_{n'l'j'}(r) r^{2+\kappa} \; dr
 * \f]
 * The part \f$ r^{2+\kappa} \f$ varies with rescaling the domain and is this
 * given by the power_kernel function.
 *
 * \tparam T        Method for calculating the wavefunction
 * \param[in] qd1   Quantum  defect data for first atom
 * \param[in] power Exponent kappa
 * \param[in] qd2   Quantum  defect data for second atom
 * \returns Radial matrix element
 */
template < typename T >
real_t IntegrateRadialElement( QuantumDefect const& qd1, int power, QuantumDefect const& qd2)
{
    T N1(qd1);
    T N2(qd2);

    auto const& x1 = N1.axis();
    auto const& y1 = N1.integrate();
    auto const& x2 = N2.axis();
    auto const& y2 = N2.integrate();
    auto const  dx = N1.dx;

    auto const xmin = x1.front() >= x2.front() ? x1.front() : x2.front();
    auto const xmax = x1.back() <= x2.back() ? x1.back() : x2.back();

    real_t mu = 0;
    // If there is an overlap, calculate the matrix element
    if (xmin <= xmax) {
        int start1 = findidx(x1, xmin);
        int end1   = findidx(x1, xmax);
        int start2 = findidx(x2, xmin);
        int end2   = findidx(x2, xmax);

        int i1, i2;
        for (i1 = start1, i2 = start2; i1 < end1 && i2 < end2; ++i1, ++i2)
        {
            mu += y1[i1]*y2[i2] * std::pow(x1[i1], T::power_kernel(power)) * dx;
        }
        mu = std::abs(2*mu);
    }

    return mu;
}


#endif // WAVEFUNCTION_H
