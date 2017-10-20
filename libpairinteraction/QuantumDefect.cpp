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
#include "EmbeddedDatabase.h"
#include "SQLite.h"

#include <algorithm>
#include <cmath>
#include <exception>
#include <iostream>
#include <sstream>
#include <string>

struct no_defect : public std::exception {
private:
    std::string msg;

public:
    no_defect(QuantumDefect const &qd) : msg()
    {
        msg = "There is no defect available for " + qd.species +
              ", n = " + std::to_string(qd.n) +
              ", l = " + std::to_string(qd.l) + ", j = " + std::to_string(qd.j);
    }

    const char *what() const throw() { return msg.c_str(); }
};

QuantumDefect::QuantumDefect(std::string const &_species, int _n, int _l,
                             double _j, std::nullptr_t)
    : ac_(), Z_(), a1_(), a2_(), a3_(), a4_(), rc_(), nstar_(), energy_(),
      species(_species), n(_n), l(_l), j(_j), ac(ac_), Z(Z_), a1(a1_), a2(a2_),
      a3(a3_), a4(a4_), rc(rc_), nstar(nstar_), energy(energy_)
{
}

QuantumDefect::QuantumDefect(std::string const &species, int n, int l, double j)
    : QuantumDefect(species, n, l, j, nullptr)
{
    static thread_local EmbeddedDatabase embedded_database{};
    setup(embedded_database);
}

QuantumDefect::QuantumDefect(std::string const &species, int n, int l, double j,
                             std::string const &database)
    : QuantumDefect(species, n, l, j, nullptr)
{
    sqlite::handle db(database, SQLITE_OPEN_READONLY);
    setup(db);
}

void QuantumDefect::setup(sqlite3 *db)
{
    std::stringstream ss;
    sqlite::statement stmt(db);
    int pot_max_l, ryd_max_l;
    int pot_l, ryd_l;
    double ryd_max_j;
    double ryd_j;

    // Determine maximal L for model potentials
    stmt.set("select MAX(L) from model_potential where (element = ?);");
    stmt.prepare();
    stmt.bind(1, species);
    if (stmt.step())
        pot_max_l = stmt.get<int>(0);
    else
        throw no_defect(*this);

    // The l to be used is the minimum of the two below
    pot_l = std::min(l, pot_max_l);

    // Determine maximal L for Rydberg-Ritz coefficients
    stmt.reset();
    stmt.set("select MAX(L) from rydberg_ritz where (element = ?);");
    stmt.prepare();
    stmt.bind(1, species);
    if (stmt.step())
        ryd_max_l = stmt.get<int>(0);
    else
        throw no_defect(*this);

    // The l to be used is the minimum of the two below
    ryd_l = std::min(l, ryd_max_l);

    // Determine maximal J for Rydberg-Ritz coefficients
    stmt.reset();
    stmt.set("select MAX(J) from rydberg_ritz where  ((element = ?1) and (L = "
             "?2));");
    stmt.prepare();
    stmt.bind(1, species);
    stmt.bind(2, ryd_l);
    if (stmt.step())
        ryd_max_j = stmt.get<double>(0);
    else
        throw no_defect(*this);

    // The j to be used is the minimum of the two below
    ryd_j = std::min(j, ryd_max_j);

    // Load model potentials from database
    stmt.reset();
    stmt.set("select ac,Z,a1,a2,a3,a4,rc from model_potential where ((element "
             "= ?1) and (L = ?2));");
    stmt.prepare();
    stmt.bind(1, species);
    stmt.bind(2, pot_l);
    if (stmt.step()) {
        ac_ = stmt.get<double>(0);
        Z_ = stmt.get<int>(1);
        a1_ = stmt.get<double>(2);
        a2_ = stmt.get<double>(3);
        a3_ = stmt.get<double>(4);
        a4_ = stmt.get<double>(5);
        rc_ = stmt.get<double>(6);
    } else
        throw no_defect(*this);

    // Load Rydberg-Ritz coefficients from database
    stmt.reset();
    stmt.set("select d0,d2,d4,d6,d8,Ry from rydberg_ritz where ((element = ?1) "
             "and (L = ?2) and (J = ?3));");
    stmt.prepare();
    stmt.bind(1, species);
    stmt.bind(2, ryd_l);
    stmt.bind(3, ryd_j);
    nstar_ = n;
    double Ry_inf = 109737.31568525;
    double Ry = Ry_inf;
    if (stmt.step()) {
        double d0 = stmt.get<double>(0);
        double d2 = stmt.get<double>(1);
        double d4 = stmt.get<double>(2);
        double d6 = stmt.get<double>(3);
        double d8 = stmt.get<double>(4);
        Ry = stmt.get<double>(5);
        nstar_ -= d0 + d2 / pow(n - d0, 2) + d4 / pow(n - d0, 4) +
                  d6 / pow(n - d0, 6) + d8 / pow(n - d0, 8);
    } else
        throw no_defect(*this);

    energy_ = -.5 * (Ry / Ry_inf) / (nstar_ * nstar_) * au2GHz;
}

double energy_level(std::string const &species, int n, int l, double j)
{
    QuantumDefect qd(species, n, l, j);
    return qd.energy;
}
