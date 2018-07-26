/*
 * Copyright (c) 2018 Sebastian Weber, Henri Menke. All rights reserved.
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

#include "MatrixElementCache.h"
#include "QuantumDefect.h"
#include "SQLite.h"
#include "version.h"

#include <cctype>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>

bool selectionRulesMomentumNew(StateOne const &state1, StateOne const &state2, int q) {
    bool validL = state1.l == state2.l;
    bool validJ = fabs(state1.j - state2.j) <= 1;
    bool validM = state1.m == state2.m + q;
    bool validQ = abs(q) <= 1;
    return validL && validJ && validM && validQ;
}

bool selectionRulesMomentumNew(StateOne const &state1, StateOne const &state2) {
    bool validL = state1.l == state2.l;
    bool validJ = fabs(state1.j - state2.j) <= 1;
    bool validM = (fabs(state1.m - state2.m) <= 1);
    return validL && validJ && validM;
}

bool selectionRulesMultipoleNew(StateOne const &state1, StateOne const &state2, int kappa, int q) {
    bool validL =
        (abs(state1.l - state2.l) <= kappa) && (kappa % 2 == abs(state1.l - state2.l) % 2);
    bool validJ = (fabs(state1.j - state2.j) <= kappa) && (state1.j + state2.j >= kappa);
    bool validM = state1.m == state2.m + q;
    bool validQ = abs(q) <= kappa;
    bool noZero = !(kappa == 2 && state1.j == state2.j && state2.j == 1.5 &&
                    state1.m == -state2.m && fabs(state1.m - state2.m) == 1);
    return validL && validJ && validM && validQ && noZero;
}

bool selectionRulesMultipoleNew(StateOne const &state1, StateOne const &state2, int kappa) {
    bool validL =
        (abs(state1.l - state2.l) <= kappa) && (kappa % 2 == abs(state1.l - state2.l) % 2);
    bool validJ = (fabs(state1.j - state2.j) <= kappa) && (state1.j + state2.j >= kappa);
    bool validM = (fabs(state1.m - state2.m) <= kappa);
    bool noZero = !(kappa == 2 && state1.j == state2.j && state2.j == 1.5 &&
                    state1.m == -state2.m && fabs(state1.m - state2.m) == 1);
    return validL && validJ && validM && noZero;
}

////////////////////////////////////////////////////////////////////
/// Constructors ///////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////

MatrixElementCache::MatrixElementCache()
    : method(NUMEROV), defectdbname(""), dbname(""), db(dbname), stmt(db),
      pid_which_created_db(utils::get_pid()) { // db and stmt have to be initialized here since they
                                               // have no default constructor
}

MatrixElementCache::MatrixElementCache(std::string const &cachedir)
    : method(NUMEROV), defectdbname(""), dbname((boost::filesystem::absolute(cachedir) /
                                                 ("cache_elements_" + version::cache() + ".db"))
                                                    .string()),
      db(dbname), stmt(db), pid_which_created_db(utils::get_pid()) {

    // Speed up database access
    stmt.exec("PRAGMA synchronous = OFF");     // do not wait on write, hand off to OS and carry on
    stmt.exec("PRAGMA journal_mode = MEMORY"); // keep rollback journal in memory during transaction

    // Create cache tables (reduced_momentum_s and reduced_momentum_l need not to be cached since
    // they are trivial)
    stmt.exec("create table if not exists cache_radial ("
              "method int, species text, k integer, n1 integer, l1 integer, j1 double,"
              "n2 integer, l2 integer, j2 double, value double, primary key (method, species, k, "
              "n1, l1, j1, n2, l2, j2)) without rowid;");

    stmt.exec(
        "create table if not exists cache_angular ("
        "k integer, j1 double, m1 double,"
        "j2 double, m2 double, value double, primary key (k, j1, m1, j2, m2)) without rowid;");

    stmt.exec(
        "create table if not exists cache_reduced_commutes_s ("
        "s double, k integer, l1 integer, j1 double,"
        "l2 integer, j2 double, value double, primary key (s, k, l1, j1, l2, j2)) without rowid;");

    stmt.exec(
        "create table if not exists cache_reduced_commutes_l ("
        "s double, k integer, l1 integer, j1 double,"
        "l2 integer, j2 double, value double, primary key (s, k, l1, j1, l2, j2)) without rowid;");

    stmt.exec("create table if not exists cache_reduced_multipole ("
              "k integer, l1 integer,"
              "l2 integer, value double, primary key (k, l1, l2)) without rowid;");
}

void MatrixElementCache::setDefectDB(std::string const &path) {
    // Set path to user defined database for quantum defects and model potential parameters
    defectdbname = path;

    // Prevent using the cache file to avoid inconsistencies through the user defined database
    dbname = "";
}

void MatrixElementCache::setMethod(method_t const &m) { method = m; }

////////////////////////////////////////////////////////////////////
/// Keys ///////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////

MatrixElementCache::CacheKey_cache_radial::CacheKey_cache_radial(method_t method,
                                                                 std::string species, int kappa,
                                                                 int n1, int n2, int l1, int l2,
                                                                 float j1, float j2)
    : species(std::move(species)), method(method), kappa(kappa) {
    if ((n1 < n2) || ((n1 == n2) && ((l1 < l2) || ((l1 == l2) && (j1 <= j2))))) {
        n = {{n1, n2}};
        l = {{l1, l2}};
        j = {{j1, j2}};
    } else {
        n = {{n2, n1}};
        l = {{l2, l1}};
        j = {{j2, j1}};
    }
}

MatrixElementCache::CacheKey_cache_angular::CacheKey_cache_angular(int kappa, float j1, float j2,
                                                                   float m1, float m2)
    : kappa(kappa) {
    if ((j1 < j2) || ((j1 == j2) && (m1 <= m2))) {
        j = {{j1, j2}};
        m = {{m1, m2}};
        sgn = 1;
    } else {
        j = {{j2, j1}};
        m = {{m2, m1}};
        sgn = pow(-1, int(j1 - m1 + j2 - m2));
    }
}

MatrixElementCache::CacheKey_cache_reduced_commutes::CacheKey_cache_reduced_commutes(
    float s, int kappa, int l1, int l2, float j1, float j2)
    : s(s), kappa(kappa) {
    if ((l1 < l2) || ((l1 == l2) && (j1 <= j2))) {
        l = {{l1, l2}};
        j = {{j1, j2}};
        sgn = 1;
    } else {
        l = {{l2, l1}};
        j = {{j2, j1}};
        sgn = pow(-1, int(l1 + j1 + l2 + j2 + 2 * s)); // TODO is this formula always correct?
    }
}

MatrixElementCache::CacheKey_cache_reduced_multipole::CacheKey_cache_reduced_multipole(int kappa,
                                                                                       int l1,
                                                                                       int l2)
    : kappa(kappa) {
    if (l1 <= l2) {
        l = {{l1, l2}};
        sgn = 1;
    } else {
        l = {{l2, l1}};
        sgn = pow(-1, kappa);
    }
}

MatrixElementCache::CacheKey_cache_radial::CacheKey_cache_radial() =
    default; // TODO the default constructors seem to be needed for serialization, can one somehow
             // circumvent the need of default constructors?

MatrixElementCache::CacheKey_cache_angular::CacheKey_cache_angular() =
    default; // TODO the default constructors seem to be needed for serialization, can one somehow
             // circumvent the need of default constructors?

MatrixElementCache::CacheKey_cache_reduced_commutes::CacheKey_cache_reduced_commutes() =
    default; // TODO the default constructors seem to be needed for serialization, can one somehow
             // circumvent the need of default constructors?

MatrixElementCache::CacheKey_cache_reduced_multipole::CacheKey_cache_reduced_multipole() =
    default; // TODO the default constructors seem to be needed for serialization, can one somehow
             // circumvent the need of default constructors?

bool MatrixElementCache::CacheKey_cache_radial::operator==(const CacheKey_cache_radial &rhs) const {
    return (method == rhs.method) && (species == rhs.species) && (kappa == rhs.kappa) &&
        (n == rhs.n) && (l == rhs.l) && (j == rhs.j);
}

bool MatrixElementCache::CacheKey_cache_angular::
operator==(const CacheKey_cache_angular &rhs) const {
    return (kappa == rhs.kappa) && (j == rhs.j) && (m == rhs.m);
}

bool MatrixElementCache::CacheKey_cache_reduced_commutes::
operator==(const CacheKey_cache_reduced_commutes &rhs) const {
    return (s == rhs.s) && (kappa == rhs.kappa) && (l == rhs.l) && (j == rhs.j);
}

bool MatrixElementCache::CacheKey_cache_reduced_multipole::
operator==(const CacheKey_cache_reduced_multipole &rhs) const {
    return (kappa == rhs.kappa) && (l == rhs.l);
}

////////////////////////////////////////////////////////////////////
/// Hasher /////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////

std::size_t MatrixElementCache::CacheKeyHasher_cache_radial::
operator()(const CacheKey_cache_radial &c) const {
    size_t seed = 0;
    boost::hash_combine(seed, c.method);
    boost::hash_combine(seed, c.species);
    boost::hash_combine(seed, c.kappa);
    boost::hash_combine(seed, c.n);
    boost::hash_combine(seed, c.l);
    boost::hash_combine(seed, c.j);
    return seed;
}

std::size_t MatrixElementCache::CacheKeyHasher_cache_angular::
operator()(const CacheKey_cache_angular &c) const {
    size_t seed = 0;
    boost::hash_combine(seed, c.kappa);
    boost::hash_combine(seed, c.j);
    boost::hash_combine(seed, c.m);
    return seed;
}

std::size_t MatrixElementCache::CacheKeyHasher_cache_reduced_commutes::
operator()(const CacheKey_cache_reduced_commutes &c) const {
    size_t seed = 0;
    boost::hash_combine(seed, c.s);
    boost::hash_combine(seed, c.kappa);
    boost::hash_combine(seed, c.l);
    boost::hash_combine(seed, c.j);
    return seed;
}

std::size_t MatrixElementCache::CacheKeyHasher_cache_reduced_multipole::
operator()(const CacheKey_cache_reduced_multipole &c) const {
    size_t seed = 0;
    boost::hash_combine(seed, c.kappa);
    boost::hash_combine(seed, c.l);
    return seed;
}

////////////////////////////////////////////////////////////////////
/// Get matrix elements ////////////////////////////////////////////
////////////////////////////////////////////////////////////////////

double MatrixElementCache::getElectricDipole(StateOne const &state_row, StateOne const &state_col) {
    return getElectricMultipole(state_row, state_col, 1, 1);
}

double MatrixElementCache::getElectricQuadrupole(StateOne const &state_row,
                                                 StateOne const &state_col) {
    return getElectricMultipole(state_row, state_col, 2, 2);
}

double MatrixElementCache::getElectricMultipole(StateOne const &state_row,
                                                StateOne const &state_col, int k) {
    return getElectricMultipole(state_row, state_col, k, k);
}

double MatrixElementCache::getDiamagnetism(StateOne const &state_row, StateOne const &state_col,
                                           int k) {
    return 2. / 3. * getElectricMultipole(state_row, state_col, 2, k);
}

double MatrixElementCache::getMagneticDipole(StateOne const &state_row, StateOne const &state_col) {
    if (state_row.species != state_col.species) {
        throw std::runtime_error("The species must be the same for the final and initial state.");
    }

    const std::string &species = state_row.species;
    const float &s = state_row.s;

    // Search cache for constituents of the matrix element
    auto key1 = CacheKey_cache_radial(method, species, 0, state_row.n, state_col.n, state_row.l,
                                      state_col.l, state_row.j, state_col.j);
    auto iter1 = cache_radial.find(key1);
    if (iter1 == cache_radial.end()) {
        cache_radial_missing.insert(key1);
    }

    auto key2 = CacheKey_cache_angular(1, state_row.j, state_col.j, state_row.m, state_col.m);
    auto iter2 = cache_angular.find(key2);
    if (iter2 == cache_angular.end()) {
        cache_angular_missing.insert(key2);
    }

    auto key3 =
        CacheKey_cache_reduced_commutes(s, 1, state_row.l, state_col.l, state_row.j, state_col.j);
    auto iter3 = cache_reduced_commutes_s.find(key3);
    if (iter3 == cache_reduced_commutes_s.end()) {
        cache_reduced_commutes_s_missing.insert(key3);
    }

    auto key4 =
        CacheKey_cache_reduced_commutes(s, 1, state_row.l, state_col.l, state_row.j, state_col.j);
    auto iter4 = cache_reduced_commutes_l.find(key4);
    if (iter4 == cache_reduced_commutes_l.end()) {
        cache_reduced_commutes_l_missing.insert(key4);
    }

    // Update cache by calculate missing constituents
    if (this->update() != 0) {
        if (iter1 == cache_radial.end()) {
            iter1 = cache_radial.find(key1);
        }
        if (iter2 == cache_angular.end()) {
            iter2 = cache_angular.find(key2);
        }
        if (iter3 == cache_reduced_commutes_s.end()) {
            iter3 = cache_reduced_commutes_s.find(key3);
        }
        if (iter4 == cache_reduced_commutes_l.end()) {
            iter4 = cache_reduced_commutes_l.find(key4);
        }
    }

    return -au2GHz / au2G * muB * iter1->second * key2.sgn * iter2->second *
        (gL * key3.sgn * iter3->second *
             sqrt(state_row.l * (state_row.l + 1) * (2 * state_row.l + 1)) +
         gS * key4.sgn * iter4->second * sqrt(s * (s + 1) * (2 * s + 1)));
}

double MatrixElementCache::getElectricMultipole(StateOne const &state_row,
                                                StateOne const &state_col, int kappa_radial,
                                                int kappa_angular) {
    if (state_row.species != state_col.species) {
        throw std::runtime_error("The species must be the same for the final and initial state.");
    }

    const std::string &species = state_row.species;
    const float &s = state_row.s;

    // Search cache for constituents of the matrix element
    auto key1 = CacheKey_cache_radial(method, species, kappa_radial, state_row.n, state_col.n,
                                      state_row.l, state_col.l, state_row.j, state_col.j);
    auto iter1 = cache_radial.find(key1);
    if (iter1 == cache_radial.end()) {
        cache_radial_missing.insert(key1);
    }

    auto key2 =
        CacheKey_cache_angular(kappa_angular, state_row.j, state_col.j, state_row.m, state_col.m);
    auto iter2 = cache_angular.find(key2);
    if (iter2 == cache_angular.end()) {
        cache_angular_missing.insert(key2);
    }

    auto key3 = CacheKey_cache_reduced_commutes(s, kappa_angular, state_row.l, state_col.l,
                                                state_row.j, state_col.j);
    auto iter3 = cache_reduced_commutes_s.find(key3);
    if (iter3 == cache_reduced_commutes_s.end()) {
        cache_reduced_commutes_s_missing.insert(key3);
    }

    auto key4 = CacheKey_cache_reduced_multipole(kappa_angular, state_row.l, state_col.l);
    auto iter4 = cache_reduced_multipole.find(key4);
    if (iter4 == cache_reduced_multipole.end()) {
        cache_reduced_multipole_missing.insert(key4);
    }

    // Update cache by calculate missing constituents
    if (this->update() != 0) {
        if (iter1 == cache_radial.end()) {
            iter1 = cache_radial.find(key1);
        }
        if (iter2 == cache_angular.end()) {
            iter2 = cache_angular.find(key2);
        }
        if (iter3 == cache_reduced_commutes_s.end()) {
            iter3 = cache_reduced_commutes_s.find(key3);
        }
        if (iter4 == cache_reduced_multipole.end()) {
            iter4 = cache_reduced_multipole.find(key4);
        }
    }

    return iter1->second * key2.sgn * iter2->second * key3.sgn * iter3->second * key4.sgn *
        iter4->second;
}

double MatrixElementCache::getRadial(StateOne const &state_row, StateOne const &state_col,
                                     int kappa) {
    if (state_row.species != state_col.species) {
        throw std::runtime_error("The species must be the same for the final and initial state.");
    }

    const std::string &species = state_row.species;

    // Search cache for constituents of the matrix element
    auto key1 = CacheKey_cache_radial(method, species, kappa, state_row.n, state_col.n, state_row.l,
                                      state_col.l, state_row.j, state_col.j);
    auto iter1 = cache_radial.find(key1);
    if (iter1 == cache_radial.end()) {
        cache_radial_missing.insert(key1);
    }

    // Update cache by calculate missing constituents
    if (this->update() != 0) {
        if (iter1 == cache_radial.end()) {
            iter1 = cache_radial.find(key1);
        }
    }

    return 1. / elementary_charge *
        iter1->second; // TODO change the converter in Wavefunction.cpp to std::pow(au2um, power) so
                       // that 1./elementary_charge is not necessary anymore
}

////////////////////////////////////////////////////////////////////
/// Precalculate matrix elements ///////////////////////////////////
////////////////////////////////////////////////////////////////////

// TODO improve methods for precalculation, for every method getSomething create the method
// precalculateSomething

void MatrixElementCache::precalculateMultipole(const std::vector<StateOne> &basis_one, int k) {
    int q = std::numeric_limits<int>::max();
    precalculate(basis_one, k, q, k, true, false, false);
}

void MatrixElementCache::precalculateRadial(const std::vector<StateOne> &basis_one, int k) {
    int q = std::numeric_limits<int>::max();
    precalculate(basis_one, k, q, k, false, false, true);
}

void MatrixElementCache::precalculateElectricMomentum(const std::vector<StateOne> &basis_one,
                                                      int q) {
    precalculate(basis_one, 1, q, 1, true, false, false);
}

void MatrixElementCache::precalculateMagneticMomentum(const std::vector<StateOne> &basis_one,
                                                      int q) {
    precalculate(basis_one, 1, q, 0, false, true, false);
}

void MatrixElementCache::precalculateDiamagnetism(const std::vector<StateOne> &basis_one, int k,
                                                  int q) {
    precalculate(basis_one, k, q, 2, true, false, false);
}

////////////////////////////////////////////////////////////////////
/// Utility methods ////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////

double MatrixElementCache::calcRadialElement(const QuantumDefect &qd1, int power,
                                             const QuantumDefect &qd2) {
    if (method == NUMEROV) {
        return IntegrateRadialElement<Numerov>(qd1, power, qd2);
    }
    if (method == WHITTAKER) {
        return IntegrateRadialElement<Whittaker>(qd1, power, qd2);
    }
    std::string msg("You have to provide all radial matrix elements on your own because you have "
                    "deactivated the calculation of missing radial matrix elements!");
    std::cout << msg << std::endl;
    throw std::runtime_error(msg);
}

void MatrixElementCache::precalculate(const std::vector<StateOne> &basis_one, int kappa_angular,
                                      int q, int kappa_radial, bool calcElectricMultipole,
                                      bool calcMagneticMomentum, bool calcRadial) {

    const std::string &species = basis_one[0].species;
    const float &s = basis_one[0].s;

    // --- Determine elements ---

    for (size_t idx_col = 0; idx_col < basis_one.size(); ++idx_col) {
        const auto &state_col = basis_one[idx_col];

        for (size_t idx_row = 0; idx_row <= idx_col; ++idx_row) {
            const auto &state_row = basis_one[idx_row];

            if (q != std::numeric_limits<int>::max() && state_row.m - state_col.m != q) {
                continue;
            }

            if (state_row.species.empty() || state_col.species.empty()) {
                continue; // TODO artifical states !!!
            }

            if ((calcElectricMultipole &&
                 selectionRulesMultipoleNew(state_row, state_col, kappa_angular)) ||
                (calcMagneticMomentum && selectionRulesMomentumNew(state_row, state_col)) ||
                calcRadial) {
                if (calcElectricMultipole || calcMagneticMomentum || calcRadial) {
                    auto key = CacheKey_cache_radial(method, species, kappa_radial, state_row.n,
                                                     state_col.n, state_row.l, state_col.l,
                                                     state_row.j, state_col.j);
                    if (cache_radial.find(key) == cache_radial.end()) {
                        cache_radial_missing.insert(key);
                    }
                }

                if (calcElectricMultipole || calcMagneticMomentum) {
                    auto key = CacheKey_cache_angular(kappa_angular, state_row.j, state_col.j,
                                                      state_row.m, state_col.m);
                    if (cache_angular.find(key) == cache_angular.end()) {
                        cache_angular_missing.insert(key);
                    }
                }

                if (calcElectricMultipole || calcMagneticMomentum) {
                    auto key = CacheKey_cache_reduced_commutes(
                        s, kappa_angular, state_row.l, state_col.l, state_row.j, state_col.j);
                    if (cache_reduced_commutes_s.find(key) == cache_reduced_commutes_s.end()) {
                        cache_reduced_commutes_s_missing.insert(key);
                    }
                }

                if (calcMagneticMomentum) {
                    auto key = CacheKey_cache_reduced_commutes(
                        s, kappa_angular, state_row.l, state_col.l, state_row.j, state_col.j);
                    if (cache_reduced_commutes_l.find(key) == cache_reduced_commutes_l.end()) {
                        cache_reduced_commutes_l_missing.insert(key);
                    }
                }

                if (calcElectricMultipole) {
                    auto key =
                        CacheKey_cache_reduced_multipole(kappa_angular, state_row.l, state_col.l);
                    if (cache_reduced_multipole.find(key) == cache_reduced_multipole.end()) {
                        cache_reduced_multipole_missing.insert(key);
                    }
                }
            }
        }
    }
}

int MatrixElementCache::update() {

    // --- Return if the cache is already up-to-date ---
    if (cache_radial_missing.empty() && cache_angular_missing.empty() &&
        cache_reduced_commutes_s_missing.empty() && cache_reduced_commutes_l_missing.empty() &&
        cache_reduced_multipole_missing.empty()) {
        return 0;
    }

    // --- Load from database ---
    if (!dbname.empty()) {

        // Reopen the connection to the database if it was opend by a different process (without
        // doing this, there can be problems with Python multiprocessing)
        if (pid_which_created_db != static_cast<long>(utils::get_pid())) {
            db = sqlite::handle(dbname);
            stmt = sqlite::statement(db);
            pid_which_created_db = utils::get_pid();

            // Speed up database access
            stmt.exec(
                "PRAGMA synchronous = OFF"); // do not wait on write, hand off to OS and carry on
            stmt.exec("PRAGMA journal_mode = MEMORY"); // keep rollback journal in memory during
                                                       // transaction
        }

        if (!cache_radial_missing.empty()) {
            stmt.set("select value from cache_radial where `method` = ?1 and `species` = ?2 and "
                     "`k` = ?3 and `n1` = ?4 and `l1` = ?5 and `j1` = ?6 and `n2` = ?7 and `l2` = "
                     "?8 and `j2` = ?9;");
            stmt.prepare();

            for (auto cached = cache_radial_missing.begin();
                 cached != cache_radial_missing.end();) {
                stmt.bind(1, cached->method);
                stmt.bind(2, cached->species);
                stmt.bind(3, cached->kappa);
                stmt.bind(4, cached->n[0]);
                stmt.bind(5, cached->l[0]);
                stmt.bind(6, cached->j[0]);
                stmt.bind(7, cached->n[1]);
                stmt.bind(8, cached->l[1]);
                stmt.bind(9, cached->j[1]);
                if (stmt.step()) {
                    cache_radial.insert({*cached, stmt.get<double>(0)});
                    cached = cache_radial_missing.erase(cached);
                } else {
                    ++cached;
                }
                stmt.reset();
            }
        }

        if (!cache_angular_missing.empty()) {
            stmt.set("select value from cache_angular where `k` = ?1 and `j1` = ?2 and `m1` = ?3 "
                     "and `j2` = ?4 and `m2` = ?5;");
            stmt.prepare();

            for (auto cached = cache_angular_missing.begin();
                 cached != cache_angular_missing.end();) {
                stmt.bind(1, cached->kappa);
                stmt.bind(2, cached->j[0]);
                stmt.bind(3, cached->m[0]);
                stmt.bind(4, cached->j[1]);
                stmt.bind(5, cached->m[1]);
                if (stmt.step()) {
                    cache_angular.insert({*cached, stmt.get<double>(0)});
                    cached = cache_angular_missing.erase(cached);
                } else {
                    ++cached;
                }
                stmt.reset();
            }
        }

        if (!cache_reduced_commutes_s_missing.empty()) {
            stmt.set("select value from cache_reduced_commutes_s where `s` = ?1 and `k` = ?2 and "
                     "`l1` = ?3 and `j1` = ?4 and `l2` = ?5 and `j2` = ?6;");
            stmt.prepare();

            for (auto cached = cache_reduced_commutes_s_missing.begin();
                 cached != cache_reduced_commutes_s_missing.end();) {
                stmt.bind(1, cached->s);
                stmt.bind(2, cached->kappa);
                stmt.bind(3, cached->l[0]);
                stmt.bind(4, cached->j[0]);
                stmt.bind(5, cached->l[1]);
                stmt.bind(6, cached->j[1]);
                if (stmt.step()) {
                    cache_reduced_commutes_s.insert({*cached, stmt.get<double>(0)});
                    cached = cache_reduced_commutes_s_missing.erase(cached);
                } else {
                    ++cached;
                }
                stmt.reset();
            }
        }

        if (!cache_reduced_commutes_l_missing.empty()) {
            stmt.set("select value from cache_reduced_commutes_l where `s` = ?1 and `k` = ?2 and "
                     "`l1` = ?3 and `j1` = ?4 and `l2` = ?5 and `j2` = ?6;");
            stmt.prepare();

            for (auto cached = cache_reduced_commutes_l_missing.begin();
                 cached != cache_reduced_commutes_l_missing.end();) {
                stmt.bind(1, cached->s);
                stmt.bind(2, cached->kappa);
                stmt.bind(3, cached->l[0]);
                stmt.bind(4, cached->j[0]);
                stmt.bind(5, cached->l[1]);
                stmt.bind(6, cached->j[1]);
                if (stmt.step()) {
                    cache_reduced_commutes_l.insert({*cached, stmt.get<double>(0)});
                    cached = cache_reduced_commutes_l_missing.erase(cached);
                } else {
                    ++cached;
                }
                stmt.reset();
            }
        }

        if (!cache_reduced_multipole_missing.empty()) {
            stmt.set("select value from cache_reduced_multipole where `k` = ?1 and `l1` = ?2 and "
                     "`l2` = ?3;");
            stmt.prepare();

            for (auto cached = cache_reduced_multipole_missing.begin();
                 cached != cache_reduced_multipole_missing.end();) {
                stmt.bind(1, cached->kappa);
                stmt.bind(2, cached->l[0]);
                stmt.bind(3, cached->l[1]);
                if (stmt.step()) {
                    cache_reduced_multipole.insert({*cached, stmt.get<double>(0)});
                    cached = cache_reduced_multipole_missing.erase(cached);
                } else {
                    ++cached;
                }
                stmt.reset();
            }
        }
    }

    // --- Calculate missing elements and write them to the database --- // TODO parallelize

    if (!dbname.empty()) {
        stmt.exec("begin transaction;");
    }

    if (!cache_radial_missing.empty()) {
        if (!dbname.empty()) {
            stmt.set("insert or ignore into cache_radial (method, species, k, n1, l1, j1, n2, l2, "
                     "j2, value) values (?1, ?2, ?3, ?4, ?5, ?6, ?7, ?8, ?9, ?10);");
            stmt.prepare();
        }

        for (auto &cached : cache_radial_missing) {
            double val = 0;
            if (defectdbname.empty()) {
                QuantumDefect qd1(cached.species, cached.n[0], cached.l[0], cached.j[0]);
                QuantumDefect qd2(cached.species, cached.n[1], cached.l[1], cached.j[1]);
                val = calcRadialElement(qd1, cached.kappa, qd2);
            } else {
                QuantumDefect qd1(cached.species, cached.n[0], cached.l[0], cached.j[0],
                                  defectdbname);
                QuantumDefect qd2(cached.species, cached.n[1], cached.l[1], cached.j[1],
                                  defectdbname);
                val = calcRadialElement(qd1, cached.kappa, qd2);
            }

            cache_radial.insert({cached, val});

            if (!dbname.empty()) {
                stmt.bind(1, cached.method);
                stmt.bind(2, cached.species);
                stmt.bind(3, cached.kappa);
                stmt.bind(4, cached.n[0]);
                stmt.bind(5, cached.l[0]);
                stmt.bind(6, cached.j[0]);
                stmt.bind(7, cached.n[1]);
                stmt.bind(8, cached.l[1]);
                stmt.bind(9, cached.j[1]);
                stmt.bind(10, val);
                stmt.step();
                stmt.reset();
            }
        }

        cache_radial_missing.clear();
    }

    if (!cache_angular_missing.empty()) {
        if (!dbname.empty()) {
            stmt.set("insert or ignore into cache_angular (k, j1, m1, j2, m2, value) values (?1, "
                     "?2, ?3, ?4, ?5, ?6);");
            stmt.prepare();
        }

        for (auto &cached : cache_angular_missing) {
            float q = cached.m[0] - cached.m[1];

            double val = pow(-1, int(cached.j[0] - cached.m[0])) *
                WignerSymbols::wigner3j(cached.j[0], cached.kappa, cached.j[1], -cached.m[0], q,
                                        cached.m[1]);
            cache_angular.insert({cached, val});

            if (!dbname.empty()) {
                stmt.bind(1, cached.kappa);
                stmt.bind(2, cached.j[0]);
                stmt.bind(3, cached.m[0]);
                stmt.bind(4, cached.j[1]);
                stmt.bind(5, cached.m[1]);
                stmt.bind(6, val);
                stmt.step();
                stmt.reset();
            }
        }

        cache_angular_missing.clear();
    }

    if (!cache_reduced_commutes_s_missing.empty()) {
        if (!dbname.empty()) {
            stmt.set("insert or ignore into cache_reduced_commutes_s (s, k, l1, j1, l2, j2, value) "
                     "values (?1, ?2, ?3, ?4, ?5, ?6, ?7);");
            stmt.prepare();
        }

        for (auto &cached : cache_reduced_commutes_s_missing) {

            // Check triangle conditions (violated e.g. for WignerSymbols::wigner6j(0, 0.5, 0.5,
            // 0.5, 0, 1))
            double val = 0;
            if (cached.l[0] >= std::abs(cached.j[0] - cached.s) &&
                cached.l[0] <= cached.j[0] + cached.s &&
                cached.l[0] >= std::abs(cached.l[1] - cached.kappa) &&
                cached.l[0] <= cached.l[1] + cached.kappa &&
                cached.j[1] >= std::abs(cached.j[0] - cached.kappa) &&
                cached.j[1] <= cached.j[0] + cached.kappa &&
                cached.j[1] >= std::abs(cached.l[1] - cached.s) &&
                cached.j[1] <=
                    cached.l[1] + cached.s) { // TODO implement this in the wignerSymbols library
                val = pow(-1, int(cached.l[0] + cached.s + cached.j[1] + cached.kappa)) *
                    sqrt((2 * cached.j[0] + 1) * (2 * cached.j[1] + 1)) *
                    WignerSymbols::wigner6j(cached.l[0], cached.j[0], cached.s, cached.j[1],
                                            cached.l[1], cached.kappa);
            }
            cache_reduced_commutes_s.insert({cached, val});

            if (!dbname.empty()) {
                stmt.bind(1, cached.s);
                stmt.bind(2, cached.kappa);
                stmt.bind(3, cached.l[0]);
                stmt.bind(4, cached.j[0]);
                stmt.bind(5, cached.l[1]);
                stmt.bind(6, cached.j[1]);
                stmt.bind(7, val);
                stmt.step();
                stmt.reset();
            }
        }

        cache_reduced_commutes_s_missing.clear();
    }

    if (!cache_reduced_commutes_l_missing.empty()) {
        if (!dbname.empty()) {
            stmt.set("insert or ignore into cache_reduced_commutes_l (s, k, l1, j1, l2, j2, value) "
                     "values (?1, ?2, ?3, ?4, ?5, ?6, ?7);");
            stmt.prepare();
        }

        for (auto &cached : cache_reduced_commutes_l_missing) {

            // Check triangle conditions (violated e.g. for WignerSymbols::wigner6j(0, 0.5, 0.5,
            // 0.5, 0, 1))
            double val = 0;
            if (cached.s >= std::abs(cached.j[0] - cached.l[0]) &&
                cached.s <= cached.j[0] + cached.l[0] &&
                cached.s >= std::abs(cached.s - cached.kappa) &&
                cached.s <= cached.s + cached.kappa &&
                cached.j[1] >= std::abs(cached.j[0] - cached.kappa) &&
                cached.j[1] <= cached.j[0] + cached.kappa &&
                cached.j[1] >= std::abs(cached.s - cached.l[0]) &&
                cached.j[1] <=
                    cached.s + cached.l[0]) { // TODO implement this in the wignerSymbols library
                val = pow(-1, int(cached.l[0] + cached.s + cached.j[0] + cached.kappa)) *
                    sqrt((2 * cached.j[0] + 1) * (2 * cached.j[1] + 1)) *
                    WignerSymbols::wigner6j(cached.s, cached.j[0], cached.l[0], cached.j[1],
                                            cached.s, cached.kappa);
            }
            cache_reduced_commutes_l.insert({cached, val});

            if (!dbname.empty()) {
                stmt.bind(1, cached.s);
                stmt.bind(2, cached.kappa);
                stmt.bind(3, cached.l[0]);
                stmt.bind(4, cached.j[0]);
                stmt.bind(5, cached.l[1]);
                stmt.bind(6, cached.j[1]);
                stmt.bind(7, val);
                stmt.step();
                stmt.reset();
            }
        }

        cache_reduced_commutes_l_missing.clear();
    }

    if (!cache_reduced_multipole_missing.empty()) {
        if (!dbname.empty()) {
            stmt.set("insert or ignore into cache_reduced_multipole (k, l1, l2, value) values (?1, "
                     "?2, ?3, ?4);");
            stmt.prepare();
        }

        for (auto &cached : cache_reduced_multipole_missing) {

            double val = pow(-1, cached.l[0]) *
                sqrt((2 * cached.l[0] + 1) * (2 * cached.l[1] + 1)) *
                WignerSymbols::wigner3j(
                             cached.l[0], cached.kappa, cached.l[1], 0, 0,
                             0); // TODO call WignerSymbols::wigner3j(cached.kappa, cached.l[1], 0,
                                 // 0, 0) and loop over the resulting vector
            cache_reduced_multipole.insert({cached, val});

            if (!dbname.empty()) {
                stmt.bind(1, cached.kappa);
                stmt.bind(2, cached.l[0]);
                stmt.bind(3, cached.l[1]);
                stmt.bind(4, val);
                stmt.step();
                stmt.reset();
            }
        }

        cache_reduced_multipole_missing.clear();
    }

    if (!dbname.empty()) {
        stmt.exec("commit transaction;");
    }

    return 1;
}

size_t MatrixElementCache::size() {
    return cache_radial.size() + cache_angular.size() + cache_reduced_commutes_s.size() +
        cache_reduced_commutes_l.size() + cache_reduced_multipole.size();
}
