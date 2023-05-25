/*
 * Copyright (c) 2016 Sebastian Weber, Henri Menke. All rights reserved.
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

#ifndef MATRIXELEMENTCACHE_H
#define MATRIXELEMENTCACHE_H

#include "Basisnames.hpp"
#include "State.hpp"
#include "Wavefunction.hpp"
#include "utils.hpp"

// clang-format off
#if __has_include (<boost/serialization/version.hpp>)
#    include <boost/serialization/version.hpp>
#endif
#if __has_include (<boost/serialization/library_version_type.hpp>)
#    include <boost/serialization/library_version_type.hpp>
#endif
// clang-format on
#include <boost/serialization/unordered_map.hpp>
#include <boost/serialization/unordered_set.hpp>
#include <wignerSymbols/wignerSymbols-cpp.h>

#include <memory>
#include <sstream>
#include <string>
#include <tuple>
#include <unordered_map>

enum method_t {
    NUMEROV = 0,
    WHITTAKER = 1,
};

class StateOne;
class StateTwo;

bool selectionRulesMomentumNew(
    StateOne const &state1, StateOne const &state2,
    int q); // TODO rename, integrate into the MatrixElementCache namespace
bool selectionRulesMomentumNew(StateOne const &state1, StateOne const &state2);
bool selectionRulesMultipoleNew(StateOne const &state1, StateOne const &state2, int kappa, int q);
bool selectionRulesMultipoleNew(StateOne const &state1, StateOne const &state2, int kappa);

class MatrixElementCache {
public:
    MatrixElementCache();
    MatrixElementCache(std::string const &cachedir);

    double getElectricDipole(StateOne const &state_row,
                             StateOne const &state_col); // return value in GHz/(V/cm)
    double getElectricMultipole(StateOne const &state_row, StateOne const &state_col,
                                int k); // return value in GHz/(V/cm)*um^(kappa-1)
    double getDiamagnetism(StateOne const &state_row, StateOne const &state_col,
                           int k); // return value in GHz/(V/cm)*um
    double getMagneticDipole(StateOne const &state_row,
                             StateOne const &state_col); // return value in GHz/G
    double
    getElectricMultipole(StateOne const &state_row, StateOne const &state_col, int kappa_radial,
                         int kappa_angular); // return value in GHz/(V/cm)*um^(kappa_radial-1)
    double getRadial(StateOne const &state_row, StateOne const &state_col,
                     int kappa); // return value in um^kappa

    void precalculateElectricMomentum(const std::vector<StateOne> &basis_one, int q);
    void precalculateMagneticMomentum(const std::vector<StateOne> &basis_one, int q);
    void precalculateDiamagnetism(const std::vector<StateOne> &basis_one, int k, int q);
    void precalculateMultipole(const std::vector<StateOne> &basis_one, int k);
    void precalculateRadial(const std::vector<StateOne> &basis_one, int k);

    void setDefectDB(std::string const &path);
    const std::string &getDefectDB() const;
    void setMethod(method_t const &m);
    void loadElectricDipoleDB(std::string const &path, std::string const &species);

    size_t size();

private:
    int update();
    void precalculate(std::shared_ptr<const BasisnamesOne> basis_one, int kappa, int q, int kappar,
                      bool calcMultipole, bool calcMomentum, bool calcRadial);
    double calcRadialElement(const QuantumDefect &qd1, int power, const QuantumDefect &qd2);
    void precalculate(const std::vector<StateOne> &basis_one, int kappa_angular, int q,
                      int kappa_radial, bool calcElectricMultipole, bool calcMagneticMomentum,
                      bool calcRadial);

    ////////////////////////////////////////////////////////////////////
    /// Keys ///////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////

    struct CacheKey_cache_radial {
        CacheKey_cache_radial(method_t method, std::string species, int kappa, int n1, int n2,
                              int l1, int l2, float j1, float j2);
        CacheKey_cache_radial();
        bool operator==(const CacheKey_cache_radial &rhs) const;
        std::string species;
        method_t method;
        int kappa;
        std::array<int, 2> n, l;
        std::array<float, 2> j;

    private:
        friend class boost::serialization::access;
        template <class Archive>
        void serialize(Archive &ar, const unsigned int /*version*/) {
            ar &method &species &kappa &n &l &j;
        }
    };

    struct CacheKey_cache_angular {
        CacheKey_cache_angular(int kappa, float j1, float j2, float m1, float m2);
        CacheKey_cache_angular();
        bool operator==(const CacheKey_cache_angular &rhs) const;
        int kappa;
        std::array<float, 2> j, m;
        int sgn;

    private:
        friend class boost::serialization::access;
        template <class Archive>
        void serialize(Archive &ar, const unsigned int /*version*/) {
            ar &kappa &j &m &sgn;
        }
    };

    struct CacheKey_cache_reduced_commutes {
        CacheKey_cache_reduced_commutes(float s, int kappa, int l1, int l2, float j1, float j2);
        CacheKey_cache_reduced_commutes();
        bool operator==(const CacheKey_cache_reduced_commutes &rhs) const;
        float s;
        int kappa;
        std::array<int, 2> l;
        std::array<float, 2> j;
        int sgn;

    private:
        friend class boost::serialization::access;
        template <class Archive>
        void serialize(Archive &ar, const unsigned int /*version*/) {
            ar &s &kappa &l &j &sgn;
        }
    };

    struct CacheKey_cache_reduced_multipole {
        CacheKey_cache_reduced_multipole(int kappa, int l1, int l2);
        CacheKey_cache_reduced_multipole();
        bool operator==(const CacheKey_cache_reduced_multipole &rhs) const;
        int kappa;
        std::array<int, 2> l;
        int sgn;

    private:
        friend class boost::serialization::access;
        template <class Archive>
        void serialize(Archive &ar, const unsigned int /*version*/) {
            ar &kappa &l &sgn;
        }
    };

    struct CacheKeyHasher_cache_radial {
        std::size_t operator()(const CacheKey_cache_radial &c) const;
    };
    struct CacheKeyHasher_cache_angular {
        std::size_t operator()(const CacheKey_cache_angular &c) const;
    };
    struct CacheKeyHasher_cache_reduced_commutes {
        std::size_t operator()(const CacheKey_cache_reduced_commutes &c) const;
    };
    struct CacheKeyHasher_cache_reduced_multipole {
        std::size_t operator()(const CacheKey_cache_reduced_multipole &c) const;
    };

    std::unordered_map<CacheKey_cache_radial, double, CacheKeyHasher_cache_radial> cache_radial;
    std::unordered_map<CacheKey_cache_angular, double, CacheKeyHasher_cache_angular> cache_angular;
    std::unordered_map<CacheKey_cache_reduced_commutes, double,
                       CacheKeyHasher_cache_reduced_commutes>
        cache_reduced_commutes_s;
    std::unordered_map<CacheKey_cache_reduced_commutes, double,
                       CacheKeyHasher_cache_reduced_commutes>
        cache_reduced_commutes_l;
    std::unordered_map<CacheKey_cache_reduced_multipole, double,
                       CacheKeyHasher_cache_reduced_multipole>
        cache_reduced_multipole;

    std::unordered_set<CacheKey_cache_radial, CacheKeyHasher_cache_radial> cache_radial_missing;
    std::unordered_set<CacheKey_cache_angular, CacheKeyHasher_cache_angular> cache_angular_missing;
    std::unordered_set<CacheKey_cache_reduced_commutes, CacheKeyHasher_cache_reduced_commutes>
        cache_reduced_commutes_s_missing;
    std::unordered_set<CacheKey_cache_reduced_commutes, CacheKeyHasher_cache_reduced_commutes>
        cache_reduced_commutes_l_missing;
    std::unordered_set<CacheKey_cache_reduced_multipole, CacheKeyHasher_cache_reduced_multipole>
        cache_reduced_multipole_missing;

    method_t method{NUMEROV};
    std::string defectdbname;
    std::string dbname;
    std::unique_ptr<sqlite::handle> db;
    std::unique_ptr<sqlite::statement> stmt;
    long pid_which_created_db;

    ////////////////////////////////////////////////////////////////////
    /// Method for serialization ///////////////////////////////////////
    ////////////////////////////////////////////////////////////////////

    friend class boost::serialization::access;
    template <class Archive>
    void serialize(Archive &ar, const unsigned int /*version*/) {
        ar &method;
        ar &dbname;
        ar &cache_radial &cache_angular &cache_reduced_commutes_s &cache_reduced_commutes_l
            &cache_reduced_multipole;
        ar &cache_radial_missing &cache_angular_missing &cache_reduced_commutes_s_missing
            &cache_reduced_commutes_l_missing &cache_reduced_multipole_missing;

        if (Archive::is_loading::value && !dbname.empty()) {
            // Open database
            db = std::make_unique<sqlite::handle>(dbname);
            stmt = std::make_unique<sqlite::statement>(*db);
            pid_which_created_db = utils::get_pid();

            // Speed up database access
            stmt->exec(
                "PRAGMA synchronous = OFF"); // do not wait on write, hand off to OS and carry on
            stmt->exec("PRAGMA journal_mode = MEMORY"); // keep rollback journal in memory during
                                                        // transaction
        }
    }
};

#endif
