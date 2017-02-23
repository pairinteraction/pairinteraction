/*
 * Copyright (c) 2016 Sebastian Weber, Henri Menke. All rights reserved.
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
#include "SQLite.h"

#include <string>
#include <iostream>
#include <sstream>
#include <exception>
#include <cmath>
#include <algorithm>


struct no_defect : public std::exception
{
    const char* what () const throw ()
    {
        return "There is no defect available";
    }
};


QuantumDefect::QuantumDefect(std::string const& _species, int _n, int _l, real_t _j)
    : species(_species), n(_n), l(_l), j(_j),
      ac(ac_), Z(Z_), a1(a1_), a2(a2_), a3(a3_), a4(a4_),
      rc(rc_), nstar(nstar_), energy(energy_)
{
    std::stringstream ss;
    sqlite::handle db("calc/databases/quantum_defects.db", SQLITE_OPEN_READONLY);
    int pot_max_l, ryd_max_l;
    int pot_l, ryd_l;
    real_t ryd_max_j;
    real_t ryd_j;

    // Determine maximal L for model potentials
    ss.str(std::string());
    ss << "select MAX(L) from model_potential where (element = '" << species << "');";
    sqlite::result res1 = db.query(ss);
    if (res1.size() > 0)
        res1.first() >> pot_max_l;
    else throw no_defect();

    // The l to be used is the minimum of the two below
    pot_l = std::min(l, pot_max_l);

    // Determine maximal L for Rydberg-Ritz coefficients
    ss.str(std::string());
    ss << "select MAX(L) from rydberg_ritz where (element = '" << species << "');";
    sqlite::result res2 = db.query(ss);
    if (res2.size() > 0)
        res2.first() >> ryd_max_l;
    else throw no_defect();

    // The l to be used is the minimum of the two below
    ryd_l = std::min(l, ryd_max_l);

    // Determine maximal J for Rydberg-Ritz coefficients
    ss.str(std::string());
    ss << "select MAX(J) from rydberg_ritz where  ("
       << "(element = '" << species << "') "
       << "and (L = " << ryd_l << ") "
       << ");";
    sqlite::result res3 = db.query(ss);
    if (res3.size() > 0)
        res3.first() >> ryd_max_j;
    else throw no_defect();

    // The j to be used is the minimum of the two below
    ryd_j = std::min(j, ryd_max_j);


    // Load model potentials from database
    ss.str(std::string());
    ss << "select ac,Z,a1,a2,a3,a4,rc from model_potential where ("
       << "(element = '" << species << "') "
       << "and (L = " << pot_l << ") "
       << ");";
    sqlite::result res4 = db.query(ss);
    if (res4.size() > 0)
        res4.first() >> ac_ >> Z_ >> a1_ >> a2_ >> a3_ >> a4_ >> rc_;
    else throw no_defect();


    // Load Rydberg-Ritz coefficients from database
    ss.str(std::string());
    ss << "select d0,d2,d4,d6,d8,Ry from rydberg_ritz where ("
       << "(element = '" << species << "') "
       << "and (L = " << ryd_l << ") "
       << "and (J = " << ryd_j << ") "
       << ");";
    sqlite::result res = db.query(ss);
    nstar_ = n;
    real_t Ry_inf = 109737.31568525; // TODO kann man hier wirklich immer den selben Wert verwenden?
    real_t d0, d2, d4, d6, d8, Ry = Ry_inf;
    if (res.size() > 0)
    {
        res.first() >> d0 >> d2 >> d4 >> d6 >> d8 >> Ry;
        nstar_ -= d0 + d2/pow(n-d0,2) + d4/pow(n-d0,4)
            + d6/pow(n-d0,6) + d8/pow(n-d0,8);
    }
    else throw no_defect();

    energy_ = -.5*(Ry/Ry_inf)/(nstar_*nstar_);
}


real_t energy_level(std::string const& species, int n, int l, real_t j) {
    QuantumDefect qd(species, n, l, j);
    return qd.energy;
}
