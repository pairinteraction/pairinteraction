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

#include <cctype>
#include <string>
#include <vector>
#include <cmath>
#include <gsl/gsl_sf_hyperg.h>

// --- Numerov's method ---

namespace model_potential {


double V(QuantumDefect const& qd, double x) {
    double Z_l = 1 + (qd.Z - 1)*std::exp(-qd.a1*x)
               - x*(qd.a3 + qd.a4*x)*std::exp(-qd.a2*x);
    double V_c = -Z_l / x;
    double V_p = -qd.ac/(2.*x*x*x*x) * (1-std::exp(-std::pow(x/qd.rc,6)));
    double V_so = 0.0;
    if ( qd.l < 4 ) {
        double s = 0.5;
        if (std::isdigit(qd.species.back())) s = (std::atoi(&qd.species.back())-1)/2.; // TODO think of a better solution
        double alpha = 7.2973525664e-3;// ~1/137 fine-structure constant CODATA 2014
        V_so = alpha*alpha/(4 * x*x*x) * (qd.j*(qd.j+1) - qd.l*(qd.l+1) - s*(s+1));
    }
    return V_c + V_p + V_so;
}


double g(QuantumDefect const& qd, double x) {
    return (2.*qd.l+.5)*(2.*qd.l+1.5)/x + 8*x*(V(qd,x) - qd.energy/au2GHz);
}


} // namespace model_potential


Numerov::Numerov(QuantumDefect const& qd)
    : qd(qd), xy()
{
    // augmented classical turning point
    double xmin = qd.n*qd.n - qd.n*std::sqrt(qd.n*qd.n-(qd.l-1)*(qd.l-1));
    if ( xmin < 2.08 )
        xmin = std::floor( std::sqrt( 2.08 ) );
    else
        xmin = std::floor( std::sqrt( xmin ) );

    double const xmax = std::sqrt( 2*qd.n*(qd.n+15) );
    double const nsteps = std::ceil( (xmax - xmin)/dx );

    xy = eigen_dense_double_t::Zero(nsteps,2);

    for (int i = 0; i < nsteps; ++i) {
        xy(i,0) = (xmin + i*dx);
    }
}


eigen_dense_double_t Numerov::integrate()
{
    using model_potential::g;

    int const nsteps = xy.rows();

    // Set the initial condition
    if ( (qd.n-qd.l) % 2 == 0 )
        xy(nsteps-2,1) = -1e-10;
    else
        xy(nsteps-2,1) = 1e-10;

    // Perform the integration using Numerov's scheme
    for (int i = nsteps-3; i >= 0; --i)
    {
        double A = (2. + 5./6. * dx*dx * g(qd,xy(i+1,0)*xy(i+1,0))) * xy(i+1,1);
        double B = (1. - 1./12.* dx*dx * g(qd,xy(i+2,0)*xy(i+2,0))) * xy(i+2,1);
        double C =  1. - 1./12.* dx*dx * g(qd,xy(i,0)*xy(i,0));
        xy(i,1) = (A - B)/C;
    }

    // Normalization
    double norm = 0;
    for (int i = 0; i < nsteps; ++i)
        norm += xy(i,1)*xy(i,1) * xy(i,0)*xy(i,0) * dx;
    norm = std::sqrt(2*norm);

    if ( norm > 0.0 ) {
        for (int i = 0; i < nsteps; ++i)
            xy(i,1) /= norm;
    }

    return xy;
}


// --- Whittaker method ---

namespace whittaker_functions {


double HypergeometricU(double a, double b, double z)
{
    if (z == 0) return NAN;
    return gsl_sf_hyperg_U(a,b,z);
}


double WhittakerW(double k, double m, double z)
{
    return std::exp(-.5*z)*std::pow(z,m+.5)*HypergeometricU(m-k+.5, 1+2*m, z);
}


double RadialWFWhittaker(double r, double nu, int l)
{
    return 1/std::sqrt(nu*nu * std::tgamma(nu+l+1) * std::tgamma(nu-l)) * WhittakerW(nu, l+.5, 2*r/nu);
}


} // namespace whittaker_functions


Whittaker::Whittaker(QuantumDefect const& qd)
    : qd(qd), xy()
{
    double const xmin = 1;
    double const xmax = std::sqrt(2*qd.n*(qd.n+15));
    double const nsteps = std::ceil( (xmax - xmin)/dx );

    xy.resize(nsteps,2);

    for (int i = 0; i < nsteps; ++i) {
        xy(i,0) = (xmin + i*dx)*(xmin + i*dx);
    }
}


eigen_dense_double_t Whittaker::integrate()
{
    using whittaker_functions::RadialWFWhittaker;

    // Set the sign
    int sign;
    if ( (qd.n-qd.l) % 2 == 0 )
        sign = -1;
    else
        sign = 1;

    int const nsteps = xy.rows();

    for (int i = 0; i < nsteps; ++i)
        xy(i,1) = sign * RadialWFWhittaker(xy(i,0), qd.nstar, qd.l);

    return xy;
}
