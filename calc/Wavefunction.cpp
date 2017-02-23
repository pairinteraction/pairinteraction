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

#include "QuantumDefect.h"
#include "Wavefunction.h"

#include <string>
#include <vector>
#include <cmath>
#include <gsl/gsl_sf_hyperg.h>

// --- Numerov's method ---

namespace model_potential {


real_t V(QuantumDefect const& qd, real_t x) {
    real_t Z_l = 1 + (qd.Z - 1)*std::exp(-qd.a1*x)
               - x*(qd.a3 + qd.a4*x)*std::exp(-qd.a2*x);
    real_t V_c = -Z_l / x;
    real_t V_p = -qd.ac/(2.*x*x*x*x) * (1-std::exp(-std::pow(x/qd.rc,6)));
    real_t V_so = 0.0;
    if ( qd.l < 4 ) {
        real_t alpha = 7.2973525664e-3;// ~1/137 fine-structure constant CODATA 2014
        V_so = alpha*alpha/(4 * x*x*x) * (qd.j*(qd.j+1) - qd.l*(qd.l+1) - 0.5*(1+0.5));
    }
    return V_c + V_p + V_so;
}


real_t g(QuantumDefect const& qd, real_t x) {
    return (2.*qd.l+.5)*(2.*qd.l+1.5)/x + 8*x*(V(qd,x) - qd.energy);
}


} // namespace model_potential


Numerov::Numerov(QuantumDefect const& qd)
    : qd(qd), x(), dx(0.01)
{
    // augmented classical turning point
    real_t xmin = qd.n*qd.n - qd.n*std::sqrt(qd.n*qd.n-(qd.l-1)*(qd.l-1));
    if ( xmin < 2.08 )
        xmin = std::floor( std::sqrt( 2.08 ) );
    else
        xmin = std::floor( std::sqrt( xmin ) );

    real_t const xmax = std::sqrt( 2*qd.n*(qd.n+15) );
    real_t const nsteps = std::ceil( (xmax - xmin)/dx );

    x.resize(nsteps);

    for (int i = 0; i < nsteps; ++i) {
        x[i] = (xmin + i*dx);
    }
}


std::vector<real_t> Numerov::axis() const
{
    return x;
}


std::vector<real_t> Numerov::integrate()
{
    using model_potential::g;

    int const nsteps = x.size();
    std::vector<real_t> y(nsteps,0.0);

    // Set the initial condition
    if ( (qd.n-qd.l) % 2 == 0 )
        y[nsteps-2] = -1e-10;
    else
        y[nsteps-2] = 1e-10;

    // Perform the integration using Numerov's scheme
    for (int i = nsteps-3; i >= 0; --i)
    {
        real_t A = (2. + 5./6. * dx*dx * g(qd,x[i+1]*x[i+1])) * y[i+1];
        real_t B = (1. - 1./12.* dx*dx * g(qd,x[i+2]*x[i+2])) * y[i+2];
        real_t C =  1. - 1./12.* dx*dx * g(qd,x[i]*x[i]);
        y[i] = (A - B)/C;
    }

    // Normalization
    real_t norm = 0;
    for (int i = 0; i < nsteps; ++i)
        norm += y[i]*y[i] * x[i]*x[i] * dx;
    norm = std::sqrt(2*norm);

    if ( norm > 0.0 ) {
        for (int i = 0; i < nsteps; ++i)
            y[i] /= norm;
    }

    return y;
}


// --- Whittaker method ---

namespace whittaker_functions {


real_t HypergeometricU(real_t a, real_t b, real_t z)
{
    if (z == 0) return NAN;
    return gsl_sf_hyperg_U(a,b,z);
}


real_t WhittakerW(real_t k, real_t m, real_t z)
{
    return std::exp(-.5*z)*std::pow(z,m+.5)*HypergeometricU(m-k+.5, 1+2*m, z);
}


real_t RadialWFWhittaker(real_t r, real_t nu, int l)
{
    return 1/std::sqrt(nu*nu * std::tgamma(nu+l+1) * std::tgamma(nu-l)) * WhittakerW(nu, l+.5, 2*r/nu);
}


} // namespace whittaker_functions


Whittaker::Whittaker(QuantumDefect const& qd)
    : qd(qd), x(), dx(0.01)
{
    real_t const xmin = 1;
    real_t const xmax = std::sqrt(2*qd.n*(qd.n+15));
    real_t const nsteps = std::ceil( (xmax - xmin)/dx );

    x.resize(nsteps);

    for (int i = 0; i < nsteps; ++i) {
        x[i] = (xmin + i*dx)*(xmin + i*dx);
    }
}


std::vector<real_t> Whittaker::axis() const
{
    return x;
}


std::vector<real_t> Whittaker::integrate()
{
    using whittaker_functions::RadialWFWhittaker;

    // Set the sign
    int sign;
    if ( (qd.n-qd.l) % 2 == 0 )
        sign = -1;
    else
        sign = 1;

    int const nsteps = x.size();
    std::vector<real_t> y(nsteps);

    for (int i = 0; i < nsteps; ++i)
        y[i] = sign * RadialWFWhittaker(x[i], qd.nstar, qd.l);

    return y;
}
