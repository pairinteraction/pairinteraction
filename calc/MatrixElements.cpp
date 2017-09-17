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

#include "Communication.h"
#include "MatrixElements.h"
#include "QuantumDefect.h"
#include "SQLite.h"
#include <sstream>
#include <iostream>
#include <string>
#include <limits>
#include <stdexcept>

bool selectionRulesMomentum(StateOne const& state1, StateOne const& state2, int q) {
    bool validL = state1.l == state2.l;
    bool validJ = fabs(state1.j-state2.j) <= 1;
    bool validM = state1.m == state2.m+q;
    bool validQ = abs(q) <= 1;
    return validL && validJ && validM && validQ;
}

bool selectionRulesMomentum(StateOne const& state1, StateOne const& state2) {
    bool validL = state1.l == state2.l;
    bool validJ = fabs(state1.j-state2.j) <= 1;
    bool validM = (fabs(state1.m - state2.m) <= 1);
    return validL && validJ && validM;
}

bool selectionRulesMultipole(StateOne const& state1, StateOne const& state2, int kappa, int q) {
    bool validL = (abs(state1.l-state2.l) <= kappa) && (kappa%2 == abs(state1.l-state2.l)%2);
    bool validJ = (fabs(state1.j-state2.j) <= kappa)  && (state1.j+state2.j >= kappa);
    bool validM = state1.m == state2.m+q;
    bool validQ = abs(q) <= kappa;
    bool noZero = !(kappa == 2 && state1.j == state2.j && state2.j == 1.5 && state1.m == -state2.m && fabs(state1.m - state2.m) == 1);
    return validL && validJ && validM && validQ && noZero;
}

bool selectionRulesMultipole(StateOne const& state1, StateOne const& state2, int kappa) {
    bool validL = (abs(state1.l-state2.l) <= kappa) && (kappa%2 == abs(state1.l-state2.l)%2);
    bool validJ = (fabs(state1.j-state2.j) <= kappa)  && (state1.j+state2.j >= kappa);
    bool validM = (fabs(state1.m - state2.m) <= kappa);
    bool noZero = !(kappa == 2 && state1.j == state2.j && state2.j == 1.5 && state1.m == -state2.m && fabs(state1.m - state2.m) == 1);
    return validL && validJ && validM && noZero;
}


MatrixElements::MatrixElements(std::string const& species, std::string const& dbname) : species(species), dbname(dbname) {
    method = "Modelpotentials";
    muB = 0.5;
    gS = 2.0023192;
    gL = 1;
}

MatrixElements::MatrixElements(const Configuration& config, std::string const& species, std::string const& dbname) : MatrixElements(species,dbname) {
    if (config["missingCalc"].str() == "true") {
        method = "Modelpotentials";
    } else if (config["missingWhittaker"].str() == "true") {
        method = "Whittaker";
    } else {
        method = "Error";
    }
}

void MatrixElements::precalculateMultipole(std::shared_ptr<const BasisnamesOne> basis_one, int k) {
    int q = std::numeric_limits<int>::max();
    precalculate(basis_one, k, q, k, true, false, false);
}

void MatrixElements::precalculateRadial(std::shared_ptr<const BasisnamesOne> basis_one, int k) {
    int q = std::numeric_limits<int>::max();
    precalculate(basis_one, k, q, k, false, false, true);
}

void MatrixElements::precalculateElectricMomentum(std::shared_ptr<const BasisnamesOne> basis_one, int q) {
    precalculate(basis_one, 1, q, 1, true, false, false);
}

void MatrixElements::precalculateMagneticMomentum(std::shared_ptr<const BasisnamesOne> basis_one, int q) {
    precalculate(basis_one, 1, q, 0, false, true, false);
}

void MatrixElements::precalculateDiamagnetism(std::shared_ptr<const BasisnamesOne> basis_one, int k, int q) {
    precalculate(basis_one, k, q, 2, true, false, false);
}

double MatrixElements::getElectricMomentum(StateOne const& state_row, StateOne const& state_col) { // TODO remove this method and use directly getMultipole(state_row,state_col,1)
    return getMultipole(state_row,state_col,1);
}

double MatrixElements::getMagneticMomentum(StateOne const& state_row, StateOne const& state_col) {
    double val = au2GHz / au2G * muB * cache_radial[0][StateTwo({{state_row.n, state_col.n}}, {{state_row.l, state_col.l}}, {{state_row.j, state_col.j}}, {{0,0}}).order()] *
            cache_angular[1][StateTwo({{0, 0}}, {{0, 0}}, {{state_row.j, state_col.j}}, {{state_row.m, state_col.m}})] *
            (gL * cache_reduced_commutes_s[1][StateTwo({{0, 0}}, {{state_row.l, state_col.l}}, {{state_row.j, state_col.j}}, {{0, 0}})] *
            sqrt(state_row.l*(state_row.l+1)*(2*state_row.l+1))  +
            gS * cache_reduced_commutes_l[1][StateTwo({{0, 0}}, {{state_row.l, state_col.l}}, {{state_row.j, state_col.j}}, {{0, 0}})] *
            sqrt(0.5*(0.5+1)*(2*0.5+1))); // TODO give the radial matrix element the unit au2GHz / au2G * muB (for the electric momentum, the radial matrix element already carries the unit)
    return val;
}

double MatrixElements::getDiamagnetism(StateOne const& state_row, StateOne const& state_col, int k) {
    double val = inverse_electron_rest_mass * 1./12.*cache_radial[2][StateTwo({{state_row.n, state_col.n}}, {{state_row.l, state_col.l}}, {{state_row.j, state_col.j}}, {{0,0}}).order()] *
            cache_angular[k][StateTwo({{0, 0}}, {{0, 0}}, {{state_row.j, state_col.j}}, {{state_row.m, state_col.m}})] *
            cache_reduced_commutes_s[k][StateTwo({{0, 0}}, {{state_row.l, state_col.l}}, {{state_row.j, state_col.j}}, {{0, 0}})] *
            cache_reduced_multipole[k][StateTwo({{0, 0}}, {{state_row.l, state_col.l}}, {{0, 0}}, {{0, 0}})]; // TODO do the multiplication with inverse_electron_rest_mass outside this class // TODO remove this method and use directly getMultipole(state_row,state_col,kappa_angular = ..., kappa_radial = 2)
    return val;
}

double MatrixElements::getMultipole(StateOne const& state_row, StateOne const& state_col, int k) {
    double val =  cache_radial[k][StateTwo({{state_row.n, state_col.n}}, {{state_row.l, state_col.l}}, {{state_row.j, state_col.j}}, {{0,0}}).order()] *
            cache_angular[k][StateTwo({{0, 0}}, {{0, 0}}, {{state_row.j, state_col.j}}, {{state_row.m, state_col.m}})] *
            cache_reduced_commutes_s[k][StateTwo({{0, 0}}, {{state_row.l, state_col.l}}, {{state_row.j, state_col.j}}, {{0, 0}})] *
            cache_reduced_multipole[k][StateTwo({{0, 0}}, {{state_row.l, state_col.l}}, {{0, 0}}, {{0, 0}})];
    return val;
}

double MatrixElements::getRadial(StateOne const& state_row, StateOne const& state_col, int k) {
    double val =  cache_radial[k][StateTwo({{state_row.n, state_col.n}}, {{state_row.l, state_col.l}}, {{state_row.j, state_col.j}}, {{0,0}}).order()];
    return val;
}

void MatrixElements::precalculate(std::shared_ptr<const BasisnamesOne> basis_one, int kappa, int q, int kappar, bool calcElectricMultipole, bool calcMagneticMomentum, bool calcRadial) {
    sqlite::handle db(dbname);

    // --- create cache tables if necessary (reduced_moemntumS and reduced_moemntumL need not to be cached since they are trivial) ---

    if (calcElectricMultipole || calcMagneticMomentum || calcRadial) {
        db.exec("CREATE TABLE IF NOT EXISTS cache_radial ("
                "method text, species text, k integer, n1 integer, l1 integer, j1 double,"
                "n2 integer, l2 integer, j2 double, value double, UNIQUE (method, species, k, n1, l1, j1, n2, l2, j2));");
        db.exec("CREATE TEMPORARY TABLE tmp_radial ("
                "n1 integer, l1 integer, j1 double,"
                "n2 integer, l2 integer, j2 double);");
    }

    if (calcElectricMultipole || calcMagneticMomentum) {
        db.exec("CREATE TABLE IF NOT EXISTS cache_angular ("
                "k integer, j1 double, m1 double,"
                "j2 double, m2 double, value double, UNIQUE (k, j1, m1, j2, m2));");
        db.exec("CREATE TEMPORARY TABLE tmp_angular ("
                "j1 double, m1 double,"
                "j2 double, m2 double);");
    }

    if (calcElectricMultipole || calcMagneticMomentum) {
        db.exec("CREATE TABLE IF NOT EXISTS cache_reduced_commutes_s ("
                "k integer, l1 integer, j1 double,"
                "l2 integer, j2 double, value double, UNIQUE (k, l1, j1, l2, j2));");
        db.exec("CREATE TEMPORARY TABLE tmp_reduced_commutes_s ("
                "l1 integer, j1 double,"
                "l2 integer, j2 double);");
    }

    if (calcMagneticMomentum) {
        db.exec("CREATE TABLE IF NOT EXISTS cache_reduced_commutes_l ("
                "k integer, l1 integer, j1 double,"
                "l2 integer, j2 double, value double, UNIQUE (k, l1, j1, l2, j2));");
        db.exec("CREATE TEMPORARY TABLE tmp_reduced_commutes_l ("
                "l1 integer, j1 double,"
                "l2 integer, j2 double);");
    }

    if (calcElectricMultipole) {
        db.exec("CREATE TABLE IF NOT EXISTS cache_reduced_multipole ("
                "k integer, l1 integer,"
                "l2 integer, value double, UNIQUE (k, l1, l2));");
        db.exec("CREATE TEMPORARY TABLE tmp_reduced_multipole ("
                "l1 integer,"
                "l2 integer);");
    }

    // --- determine elements ---

    std::stringstream ss;

    db.exec("begin transaction;");

    for (const auto &state_col : *basis_one) {
        for (const auto &state_row : *basis_one) {
            if (q != std::numeric_limits<int>::max() && state_row.m-state_col.m != q) {
                continue;
            }

            //if (state_row.idx < state_col.idx) { // TODO
            //    continue;
            //}

            if ((calcElectricMultipole && selectionRulesMultipole(state_row, state_col, kappa)) || (calcMagneticMomentum && selectionRulesMomentum(state_row, state_col)) || calcRadial) {
                if (calcElectricMultipole || calcMagneticMomentum  || calcRadial) {
                    StateTwo state_nlj = StateTwo({{state_row.n, state_col.n}}, {{state_row.l, state_col.l}}, {{state_row.j, state_col.j}}, {{0,0}}).order();
                    auto missing_cache_radial = cache_radial[kappar].insert({state_nlj, std::numeric_limits<double>::max()});
                    if (missing_cache_radial.second) {
                        ss.str(std::string());
                        ss << "insert into tmp_radial (n1,l1,j1,n2,l2,j2) values ("
                           << state_nlj.n[0] << "," << state_nlj.l[0] << "," << state_nlj.j[0] << ","
                           << state_nlj.n[1] << "," << state_nlj.l[1] << "," << state_nlj.j[1] << ");";
                        db.exec(ss.str());
                    }
                }

                if (calcElectricMultipole || calcMagneticMomentum) {
                    StateTwo state_jm = StateTwo({{0, 0}}, {{0, 0}}, {{state_row.j, state_col.j}}, {{state_row.m, state_col.m}});
                    auto missing_cache_angular = cache_angular[kappa].insert({state_jm, std::numeric_limits<double>::max()});
                    if (missing_cache_angular.second) {
                        ss.str(std::string());
                        ss << "insert into tmp_angular (j1,m1,j2,m2) values ("
                           << state_jm.j[0] << "," << state_jm.m[0] << ","
                           << state_jm.j[1] << "," << state_jm.m[1] << ");";
                        db.exec(ss.str());
                    }
                }

                if (calcElectricMultipole || calcMagneticMomentum) {
                    StateTwo state_lj = StateTwo({{0, 0}}, {{state_row.l, state_col.l}}, {{state_row.j, state_col.j}}, {{0, 0}});
                    auto missing_cache_reduced_commutes_s = cache_reduced_commutes_s[kappa].insert({state_lj, std::numeric_limits<double>::max()});
                    if (missing_cache_reduced_commutes_s.second) {
                        ss.str(std::string());
                        ss << "insert into tmp_reduced_commutes_s (l1,j1,l2,j2) values ("
                           << state_lj.l[0] << "," << state_lj.j[0] << ","
                           << state_lj.l[1] << "," << state_lj.j[1] << ");";
                        db.exec(ss.str());
                    }
                }

                if (calcMagneticMomentum) {
                    StateTwo state_lj = StateTwo({{0, 0}}, {{state_row.l, state_col.l}}, {{state_row.j, state_col.j}}, {{0, 0}});
                    auto missing_cache_reduced_commutes_l = cache_reduced_commutes_l[kappa].insert({state_lj, std::numeric_limits<double>::max()});
                    if (missing_cache_reduced_commutes_l.second) {
                        ss.str(std::string());
                        ss << "insert into tmp_reduced_commutes_l (l1,j1,l2,j2) values ("
                           << state_lj.l[0] << "," << state_lj.j[0] << ","
                           << state_lj.l[1] << "," << state_lj.j[1] << ");";
                        db.exec(ss.str());
                    }
                }

                if (calcElectricMultipole) {
                    StateTwo state_l = StateTwo({{0, 0}}, {{state_row.l, state_col.l}}, {{0, 0}}, {{0, 0}});
                    auto missing_cache_reduced_multipole = cache_reduced_multipole[kappa].insert({state_l, std::numeric_limits<double>::max()});
                    if (missing_cache_reduced_multipole.second) {
                        ss.str(std::string());
                        ss << "insert into tmp_reduced_multipole (l1,l2) values ("
                           << state_l.l[0] << ","
                           << state_l.l[1] << ");";
                        db.exec(ss.str());
                    }
                }
            }
        }
    }

    db.exec("end transaction;");

    // --- load from database ---

    int n1, n2, l1, l2;
    float j1, j2, m1, m2;
    double value;


    if (calcElectricMultipole || calcMagneticMomentum  || calcRadial) {
        ss.str(std::string());
        ss << "SELECT c.n1, c.l1, c.j1, c.n2, c.l2, c.j2, c.value FROM cache_radial c INNER JOIN tmp_radial t ON ("
           << "c.n1 = t.n1 AND c.l1 = t.l1 AND c.j1 = t.j1 AND c.n2 = t.n2 AND c.l2 = t.l2 AND c.j2 = t.j2) "
           << "WHERE c.method= '" << method << "' AND c.species = '" << species << "' AND c.k = " << kappar << ";";
        sqlite::statement stmt(db, ss.str());
        stmt.prepare();
        while (stmt.step())
        {
            n1 = stmt.get<decltype(n1)>(0);
            l1 = stmt.get<decltype(l1)>(1);
            j1 = stmt.get<decltype(j1)>(2);
            n2 = stmt.get<decltype(n2)>(3);
            l2 = stmt.get<decltype(l2)>(4);
            j2 = stmt.get<decltype(j2)>(5);
            value = stmt.get<decltype(value)>(6);
            cache_radial[kappar][StateTwo({{n1, n2}}, {{l1, l2}}, {{j1, j2}}, {{0,0}})] = value;
        }
    }

    if (calcElectricMultipole || calcMagneticMomentum) {
        ss.str(std::string());
        ss << "SELECT c.j1, c.m1, c.j2, c.m2, c.value FROM cache_angular c INNER JOIN tmp_angular t ON ("
           << "c.j1 = t.j1 AND c.m1 = t.m1 AND c.j2 = t.j2 AND c.m2 = t.m2) "
           << "WHERE c.k = " << kappa << ";";
        sqlite::statement stmt(db, ss.str());
        stmt.prepare();
        while (stmt.step())
        {
            j1 = stmt.get<decltype(j1)>(0);
            m1 = stmt.get<decltype(m1)>(1);
            j2 = stmt.get<decltype(j2)>(2);
            m2 = stmt.get<decltype(m2)>(3);
            value = stmt.get<decltype(value)>(4);
            cache_angular[kappa][StateTwo({{0, 0}}, {{0, 0}}, {{j1, j2}}, {{m1,m2}})] = value;
        }
    }

    if (calcElectricMultipole || calcMagneticMomentum) {
        ss.str(std::string());
        ss << "SELECT c.l1, c.j1, c.l2, c.j2, c.value FROM cache_reduced_commutes_s c INNER JOIN tmp_reduced_commutes_s t ON ("
           << "c.l1 = t.l1 AND c.j1 = t.j1 AND c.l2 = t.l2 AND c.j2 = t.j2) "
           << "WHERE c.k = " << kappa << ";";
        sqlite::statement stmt(db, ss.str());
        stmt.prepare();
        while (stmt.step())
        {
            l1 = stmt.get<decltype(l1)>(0);
            j1 = stmt.get<decltype(j1)>(1);
            l2 = stmt.get<decltype(l2)>(2);
            j2 = stmt.get<decltype(j2)>(3);
            value = stmt.get<decltype(value)>(4);
            cache_reduced_commutes_s[kappa][StateTwo({{0, 0}}, {{l1, l2}}, {{j1, j2}}, {{0,0}})] = value;
        }
    }

    if (calcMagneticMomentum) {
        ss.str(std::string());
        ss << "SELECT c.l1, c.j1, c.l2, c.j2, c.value FROM cache_reduced_commutes_l c INNER JOIN tmp_reduced_commutes_l t ON ("
           << "c.l1 = t.l1 AND c.j1 = t.j1 AND c.l2 = t.l2 AND c.j2 = t.j2) "
           << "WHERE c.k = " << kappa << ";";
        sqlite::statement stmt(db, ss.str());
        stmt.prepare();
        while (stmt.step())
        {
            l1 = stmt.get<decltype(l1)>(0);
            j1 = stmt.get<decltype(j1)>(1);
            l2 = stmt.get<decltype(l2)>(2);
            j2 = stmt.get<decltype(j2)>(3);
            value = stmt.get<decltype(value)>(4);
            cache_reduced_commutes_l[kappa][StateTwo({{0, 0}}, {{l1, l2}}, {{j1, j2}}, {{0,0}})] = value;
        }
    }

    if (calcElectricMultipole) {
        ss.str(std::string());
        ss << "SELECT c.l1, c.l2, c.value FROM cache_reduced_multipole c INNER JOIN tmp_reduced_multipole t ON ("
           << "c.l1 = t.l1 AND c.l2 = t.l2) "
           << "WHERE c.k = " << kappa << ";";
        sqlite::statement stmt(db, ss.str());
        stmt.prepare();
        while (stmt.step())
        {
            l1 = stmt.get<decltype(l1)>(0);
            l2 = stmt.get<decltype(l2)>(1);
            value = stmt.get<decltype(value)>(2);
            cache_reduced_multipole[kappa][StateTwo({{0, 0}}, {{l1, l2}}, {{0, 0}}, {{0,0}})] = value;
        }
    }

    // --- calculate missing elements and write them to the database ---

    db.exec("begin transaction;");

    if (calcElectricMultipole || calcMagneticMomentum  || calcRadial) {
        for (auto &cache : cache_radial[kappar]) {
            if (cache.second == std::numeric_limits<double>::max()) {
                auto state = cache.first;

                QuantumDefect qd1(species, state.n[0], state.l[0], state.j[0]);
                QuantumDefect qd2(species, state.n[1], state.l[1], state.j[1]);
                cache.second = calcRadialElement(qd1, kappar, qd2);

                ss.str(std::string());
                ss << "insert into cache_radial (method, species, k, n1, l1, j1, n2, l2, j2, value) values ("
                   << "'" << method << "'" << "," << "'" << species << "'" << "," << kappar << ","
                   << state.n[0] << "," << state.l[0] << "," << state.j[0] << ","
                   << state.n[1] << "," << state.l[1] << "," << state.j[1] << ","
                   << cache.second << ");";
                db.exec(ss.str());

            }
        }
    }

    if (calcElectricMultipole || calcMagneticMomentum) {
        for (auto &cache : cache_angular[kappa]) {
            if (cache.second == std::numeric_limits<double>::max()) {
                auto state = cache.first;
                float q = state.m[0]-state.m[1];

                cache.second = pow(-1, state.j[0]-state.m[0]) *
                        WignerSymbols::wigner3j(state.j[0], kappa, state.j[1], -state.m[0], q, state.m[1]);

                ss.str(std::string());
                ss << "insert into cache_angular (k, j1, m1, j2, m2, value) values ("
                   << kappa << ","
                   << state.j[0] << "," << state.m[0] << ","
                   << state.j[1] << "," << state.m[1] << ","
                   << cache.second << ");";
                db.exec(ss.str());
            }
        }
    }

    if (calcElectricMultipole || calcMagneticMomentum) {
        for (auto &cache : cache_reduced_commutes_s[kappa]) {
            if (cache.second == std::numeric_limits<double>::max()) {
                auto state = cache.first;

                cache.second = pow(-1, state.l[0] + 0.5 + state.j[1] + kappa) * sqrt((2*state.j[0]+1)*(2*state.j[1]+1)) *
                        WignerSymbols::wigner6j(state.l[0], state.j[0], 0.5, state.j[1], state.l[1], kappa);

                ss.str(std::string());
                ss << "insert into cache_reduced_commutes_s (k, l1, j1, l2, j2, value) values ("
                   << kappa << ","
                   << state.l[0] << "," << state.j[0] << ","
                   << state.l[1] << "," << state.j[1] << ","
                   << cache.second << ");";
                db.exec(ss.str());
            }
        }
    }

    if (calcMagneticMomentum) {
        for (auto &cache : cache_reduced_commutes_l[kappa]) {
            if (cache.second == std::numeric_limits<double>::max()) {
                auto state = cache.first;

                cache.second = pow(-1, state.l[0] + 0.5 + state.j[0] + kappa) * sqrt((2*state.j[0]+1)*(2*state.j[1]+1)) *
                        WignerSymbols::wigner6j(0.5, state.j[0], state.l[0], state.j[1], 0.5, kappa);

                ss.str(std::string());
                ss << "insert into cache_reduced_commutes_l (k, l1, j1, l2, j2, value) values ("
                   << kappa << ","
                   << state.l[0] << "," << state.j[0] << ","
                   << state.l[1] << "," << state.j[1] << ","
                   << cache.second << ");";
                db.exec(ss.str());
            }
        }
    }

    if (calcElectricMultipole) {
        for (auto &cache : cache_reduced_multipole[kappa]) {
            if (cache.second == std::numeric_limits<double>::max()) {
                auto state = cache.first;

                cache.second = pow(-1, state.l[0]) * sqrt((2*state.l[0]+1)*(2*state.l[1]+1)) *
                        WignerSymbols::wigner3j(state.l[0], kappa, state.l[1], 0, 0, 0);

                ss.str(std::string());
                ss << "insert into cache_reduced_multipole (k, l1, l2, value) values ("
                   << kappa << ","
                   << state.l[0] << ","
                   << state.l[1] << ","
                   << cache.second << ");";
                db.exec(ss.str());
            }
        }
    }

    db.exec("end transaction;");
}

double MatrixElements::calcRadialElement(const QuantumDefect &qd1, int power,
                                         const QuantumDefect &qd2) {
    if (method == "Modelpotentials") {
        return IntegrateRadialElement<Numerov>(qd1, power, qd2);
    } else if(method == "Whittaker") {
        return IntegrateRadialElement<Whittaker>(qd1, power, qd2);
    } else {
        std::string msg("You have to provide all radial matrix elements on your own because you have deactivated the calculation of missing radial matrix elements!");
        auto context = zmq::context();
        auto publisher = context.socket(ZMQ_PUB);
        publisher.bind(zmq::endpoint::name.c_str());
        publisher.send(">>ERR%s", msg.c_str()); // TODO make it thread save
        throw std::runtime_error(msg);
    }
}

















void MatrixElements::precalculateMultipole(const std::vector<StateOne> &basis_one, int k) {
    int q = std::numeric_limits<int>::max();
    precalculate(basis_one, k, q, k, true, false, false);
}

void MatrixElements::precalculateRadial(const std::vector<StateOne> &basis_one, int k) {
    int q = std::numeric_limits<int>::max();
    precalculate(basis_one, k, q, k, false, false, true);
}

void MatrixElements::precalculateElectricMomentum(const std::vector<StateOne> &basis_one, int q) {
    precalculate(basis_one, 1, q, 1, true, false, false);
}

void MatrixElements::precalculateMagneticMomentum(const std::vector<StateOne> &basis_one, int q) {
    precalculate(basis_one, 1, q, 0, false, true, false);
}

void MatrixElements::precalculateDiamagnetism(const std::vector<StateOne> &basis_one, int k, int q) {
    precalculate(basis_one, k, q, 2, true, false, false);
}

void MatrixElements::precalculate(const std::vector<StateOne> &basis_one, int kappa, int q, int kappar, bool calcElectricMultipole, bool calcMagneticMomentum, bool calcRadial) {
    sqlite::handle db(dbname);
    sqlite::statement stmt(db);

    // --- Speed up database access ---

    stmt.exec("PRAGMA synchronous = OFF"); // do not wait on write, hand off to OS and carry on
    stmt.exec("PRAGMA journal_mode = MEMORY"); // keep rollback journal in memory during transaction

    // --- Create cache tables if necessary (reduced_moemntum_s and reduced_moemntum_l need not to be cached since they are trivial) --- // TODO put this into the constructor of the prospective cache object

    if (calcElectricMultipole || calcMagneticMomentum || calcRadial) {
        stmt.exec("create table if not exists cache_radial ("
                 "method text, species text, k integer, n1 integer, l1 integer, j1 double,"
                 "n2 integer, l2 integer, j2 double, value double, primary key (method, species, k, n1, l1, j1, n2, l2, j2)) without rowid;");
    }

    if (calcElectricMultipole || calcMagneticMomentum) {
        stmt.exec("create table if not exists cache_angular ("
                 "k integer, j1 double, m1 double,"
                 "j2 double, m2 double, value double, primary key (k, j1, m1, j2, m2)) without rowid;");
    }

    if (calcElectricMultipole || calcMagneticMomentum) {
        stmt.exec("create table if not exists cache_reduced_commutes_s ("
                 "k integer, l1 integer, j1 double,"
                 "l2 integer, j2 double, value double, primary key (k, l1, j1, l2, j2)) without rowid;");
    }

    if (calcMagneticMomentum) {
        stmt.exec("create table if not exists cache_reduced_commutes_l ("
                 "k integer, l1 integer, j1 double,"
                 "l2 integer, j2 double, value double, primary key (k, l1, j1, l2, j2)) without rowid;");
    }

    if (calcElectricMultipole) {
        stmt.exec("create table if not exists cache_reduced_multipole ("
                 "k integer, l1 integer,"
                 "l2 integer, value double, primary key (k, l1, l2)) without rowid;");
    }

    // --- Determine elements ---

    for (const auto &state_col : basis_one) {
        for (const auto &state_row : basis_one) {
            if (q != std::numeric_limits<int>::max() && state_row.m-state_col.m != q) {
                continue;
            }

            if (state_row.element.empty() || state_col.element.empty()  ) continue; // TODO artifical states !!!

            //if (state_row.idx < state_col.idx) { // TODO
            //    continue;
            //}

            if ((calcElectricMultipole && selectionRulesMultipole(state_row, state_col, kappa)) || (calcMagneticMomentum && selectionRulesMomentum(state_row, state_col)) || calcRadial) {
                if (calcElectricMultipole || calcMagneticMomentum  || calcRadial) {
                    StateTwo state_nlj = StateTwo({{state_row.n, state_col.n}}, {{state_row.l, state_col.l}}, {{state_row.j, state_col.j}}, {{0,0}}).order();
                    cache_radial[kappar].insert({state_nlj, std::numeric_limits<double>::max()}); // TODO can this be speed up?
                }

                if (calcElectricMultipole || calcMagneticMomentum) {
                    StateTwo state_jm = StateTwo({{0, 0}}, {{0, 0}}, {{state_row.j, state_col.j}}, {{state_row.m, state_col.m}});
                    cache_angular[kappa].insert({state_jm, std::numeric_limits<double>::max()});
                }

                if (calcElectricMultipole || calcMagneticMomentum) {
                    StateTwo state_lj = StateTwo({{0, 0}}, {{state_row.l, state_col.l}}, {{state_row.j, state_col.j}}, {{0, 0}});
                    cache_reduced_commutes_s[kappa].insert({state_lj, std::numeric_limits<double>::max()});
                }

                if (calcMagneticMomentum) {
                    StateTwo state_lj = StateTwo({{0, 0}}, {{state_row.l, state_col.l}}, {{state_row.j, state_col.j}}, {{0, 0}});
                    cache_reduced_commutes_l[kappa].insert({state_lj, std::numeric_limits<double>::max()});
                }

                if (calcElectricMultipole) {
                    StateTwo state_l = StateTwo({{0, 0}}, {{state_row.l, state_col.l}}, {{0, 0}}, {{0, 0}});
                    cache_reduced_multipole[kappa].insert({state_l, std::numeric_limits<double>::max()});
                }
            }
        }
    }

    // --- Load from database ---

    if (calcElectricMultipole || calcMagneticMomentum  || calcRadial) {
        stmt.set("select value from cache_radial where `method` = ?1 and `species` = ?2 and `k` = ?3 and `n1` = ?4 and `l1` = ?5 and `j1` = ?6 and `n2` = ?7 and `l2` = ?8 and `j2` = ?9;");

        for (auto &cache : cache_radial[kappar]) {
            auto state = cache.first;
            stmt.prepare(); // @henri: can this be moved outside the loop?
            stmt.bind(1, method);
            stmt.bind(2, species);
            stmt.bind(3, kappar);
            stmt.bind(4, state.n[0]);
            stmt.bind(5, state.l[0]);
            stmt.bind(6, state.j[0]);
            stmt.bind(7, state.n[1]);
            stmt.bind(8, state.l[1]);
            stmt.bind(9, state.j[1]);
            if (stmt.step()) cache.second = stmt.get<double>(0);
            stmt.reset(); // @henri: is this needed here?
        }
    }

    if (calcElectricMultipole || calcMagneticMomentum) {
        stmt.set("select value from cache_angular where `k` = ?1 and `j1` = ?2 and `m1` = ?3 and `j2` = ?4 and `m2` = ?5;");

        for (auto &cache : cache_angular[kappa]) {
            auto state = cache.first;
            stmt.prepare(); // @henri: can this be moved outside the loop?
            stmt.bind(1, kappa);
            stmt.bind(2, state.j[0]);
            stmt.bind(3, state.m[0]);
            stmt.bind(4, state.j[1]);
            stmt.bind(5, state.m[1]);
            if (stmt.step()) cache.second = stmt.get<double>(0);
            stmt.reset(); // @henri: is this needed here?
        }
    }

    if (calcElectricMultipole || calcMagneticMomentum) {
        stmt.set("select value from cache_reduced_commutes_s where `k` = ?1 and `l1` = ?2 and `j1` = ?3 and `l2` = ?4 and `j2` = ?5;");

        for (auto &cache : cache_reduced_commutes_s[kappa]) {
            auto state = cache.first;
            stmt.prepare(); // @henri: can this be moved outside the loop?
            stmt.bind(1, kappa);
            stmt.bind(2, state.l[0]);
            stmt.bind(3, state.j[0]);
            stmt.bind(4, state.l[1]);
            stmt.bind(5, state.j[1]);
            if (stmt.step()) cache.second = stmt.get<double>(0);
            stmt.reset(); // @henri: is this needed here?
        }
    }

    if (calcMagneticMomentum) {
        stmt.set("select value from cache_reduced_commutes_l where `k` = ?1 and `l1` = ?2 and `j1` = ?3 and `l2` = ?4 and `j2` = ?5;");

        for (auto &cache : cache_reduced_commutes_l[kappa]) {
            auto state = cache.first;
            stmt.prepare(); // @henri: can this be moved outside the loop?
            stmt.bind(1, kappa);
            stmt.bind(2, state.l[0]);
            stmt.bind(3, state.j[0]);
            stmt.bind(4, state.l[1]);
            stmt.bind(5, state.j[1]);
            if (stmt.step()) cache.second = stmt.get<double>(0);
            stmt.reset(); // @henri: is this needed here?
        }
    }

    if (calcElectricMultipole) {
        stmt.set("select value from cache_reduced_multipole where `k` = ?1 and `l1` = ?2 and `l2` = ?3;");

        for (auto &cache : cache_reduced_multipole[kappa]) {
            auto state = cache.first;
            stmt.prepare(); // @henri: can this be moved outside the loop?
            stmt.bind(1, kappa);
            stmt.bind(2, state.l[0]);
            stmt.bind(3, state.l[1]);
            if (stmt.step()) cache.second = stmt.get<double>(0);
            stmt.reset(); // @henri: is this needed here?
        }
    }

    // --- Calculate missing elements and write them to the database ---

    stmt.exec("begin transaction;");

    if (calcElectricMultipole || calcMagneticMomentum  || calcRadial) {
        stmt.set("insert or ignore into cache_radial (method, species, k, n1, l1, j1, n2, l2, j2, value) values (?1, ?2, ?3, ?4, ?5, ?6, ?7, ?8, ?9, ?10);");

        for (auto &cache : cache_radial[kappar]) {
            if (cache.second == std::numeric_limits<double>::max()) {
                auto state = cache.first;

                QuantumDefect qd1(species, state.n[0], state.l[0], state.j[0]);
                QuantumDefect qd2(species, state.n[1], state.l[1], state.j[1]);
                cache.second = calcRadialElement(qd1, kappar, qd2);

                stmt.prepare(); // @henri: can this be moved outside the loop?
                stmt.bind(1, method);
                stmt.bind(2, species);
                stmt.bind(3, kappar);
                stmt.bind(4, state.n[0]);
                stmt.bind(5, state.l[0]);
                stmt.bind(6, state.j[0]);
                stmt.bind(7, state.n[1]);
                stmt.bind(8, state.l[1]);
                stmt.bind(9, state.j[1]);
                stmt.bind(10, cache.second);

                stmt.step();
                stmt.reset(); // @henri: can this be moved outside the loop?
            }
        }
    }

    if (calcElectricMultipole || calcMagneticMomentum) {
        stmt.set("insert or ignore into cache_angular (k, j1, m1, j2, m2, value) values (?1, ?2, ?3, ?4, ?5, ?6);");

        for (auto &cache : cache_angular[kappa]) {
            if (cache.second == std::numeric_limits<double>::max()) {
                auto state = cache.first;
                float q = state.m[0]-state.m[1];

                cache.second = pow(-1, state.j[0]-state.m[0]) *
                        WignerSymbols::wigner3j(state.j[0], kappa, state.j[1], -state.m[0], q, state.m[1]);

                if (cache.second > 1e9) {
                    std::cout << "Warning: Error in the calculation of the Wigner-3j-symbol detected and resolved." << std::endl;
                    cache.second = 0;
                }

                stmt.prepare(); // @henri: can this be moved outside the loop?
                stmt.bind(1, kappa);
                stmt.bind(2, state.j[0]);
                stmt.bind(3, state.m[0]);
                stmt.bind(4, state.j[1]);
                stmt.bind(5, state.m[1]);
                stmt.bind(6, cache.second);
                stmt.step();
                stmt.reset(); // @henri: can this be moved outside the loop?
            }
        }
    }

    if (calcElectricMultipole || calcMagneticMomentum) {
        stmt.set("insert or ignore into cache_reduced_commutes_s (k, l1, j1, l2, j2, value) values (?1, ?2, ?3, ?4, ?5, ?6);");

        for (auto &cache : cache_reduced_commutes_s[kappa]) {
            if (cache.second == std::numeric_limits<double>::max()) {
                auto state = cache.first;

                cache.second = pow(-1, state.l[0] + 0.5 + state.j[1] + kappa) * sqrt((2*state.j[0]+1)*(2*state.j[1]+1)) *
                        WignerSymbols::wigner6j(state.l[0], state.j[0], 0.5, state.j[1], state.l[1], kappa);

                if (cache.second > 1e9) {
                    std::cout << "Warning: Error in the calculation of the Wigner-6j-symbol detected and resolved." << std::endl; // happens e.g. for WignerSymbols::wigner6j(0, 0.5, 0.5, 0.5, 0, 1)
                    cache.second = 0;
                }

                stmt.prepare(); // @henri: can this be moved outside the loop?
                stmt.bind(1, kappa);
                stmt.bind(2, state.l[0]);
                stmt.bind(3, state.j[0]);
                stmt.bind(4, state.l[1]);
                stmt.bind(5, state.j[1]);
                stmt.bind(6, cache.second);
                stmt.step();
                stmt.reset(); // @henri: can this be moved outside the loop?
            }
        }
    }

    if (calcMagneticMomentum) {
        stmt.set("insert or ignore into cache_reduced_commutes_l (k, l1, j1, l2, j2, value) values (?1, ?2, ?3, ?4, ?5, ?6);");

        for (auto &cache : cache_reduced_commutes_l[kappa]) {
            if (cache.second == std::numeric_limits<double>::max()) {
                auto state = cache.first;

                cache.second = pow(-1, state.l[0] + 0.5 + state.j[0] + kappa) * sqrt((2*state.j[0]+1)*(2*state.j[1]+1)) *
                        WignerSymbols::wigner6j(0.5, state.j[0], state.l[0], state.j[1], 0.5, kappa);

                if (cache.second > 1e9) {
                    std::cout << "Warning: Error in the calculation of the Wigner-6j-symbol detected and resolved." << std::endl; // happens e.g. for WignerSymbols::wigner6j(0, 0.5, 0.5, 0.5, 0, 1)
                    cache.second = 0;
                }

                stmt.prepare(); // @henri: can this be moved outside the loop?
                stmt.bind(1, kappa);
                stmt.bind(2, state.l[0]);
                stmt.bind(3, state.j[0]);
                stmt.bind(4, state.l[1]);
                stmt.bind(5, state.j[1]);
                stmt.bind(6, cache.second);
                stmt.step();
                stmt.reset(); // @henri: can this be moved outside the loop?
            }
        }
    }

    if (calcElectricMultipole) {
        stmt.set("insert or ignore into cache_reduced_multipole (k, l1, l2, value) values (?1, ?2, ?3, ?4);");

        for (auto &cache : cache_reduced_multipole[kappa]) {
            if (cache.second == std::numeric_limits<double>::max()) {
                auto state = cache.first;

                cache.second = pow(-1, state.l[0]) * sqrt((2*state.l[0]+1)*(2*state.l[1]+1)) *
                        WignerSymbols::wigner3j(state.l[0], kappa, state.l[1], 0, 0, 0);

                if (cache.second > 1e9) {
                    std::cout << "Warning: Error in the calculation of the Wigner-3j-symbol detected and resolved." << std::endl;
                    cache.second = 0;
                }

                stmt.prepare(); // @henri: can this be moved outside the loop?
                stmt.bind(1, kappa);
                stmt.bind(2, state.l[0]);
                stmt.bind(3, state.l[1]);
                stmt.bind(4, cache.second);
                stmt.step();
                stmt.reset(); // @henri: can this be moved outside the loop?
            }
        }
    }

    stmt.exec("commit transaction;");
}
