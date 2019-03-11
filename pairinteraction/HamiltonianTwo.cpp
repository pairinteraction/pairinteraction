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

#include "HamiltonianTwo.h"
#include <stdexcept>
#include <utility>

#include <boost/format.hpp>

HamiltonianTwo::HamiltonianTwo(const Configuration &config, boost::filesystem::path &path_cache,
                               const std::shared_ptr<HamiltonianOne> &hamiltonian_one)
    : hamiltonian_one1(hamiltonian_one), hamiltonian_one2(hamiltonian_one),
      path_cache(path_cache) { // TODO

    samebasis = true;

    calculate(config);
}

HamiltonianTwo::HamiltonianTwo(const Configuration &config, boost::filesystem::path &path_cache,
                               std::shared_ptr<HamiltonianOne> hamiltonian_one1,
                               std::shared_ptr<HamiltonianOne> hamiltonian_one2)
    : hamiltonian_one1(std::move(hamiltonian_one1)), hamiltonian_one2(std::move(hamiltonian_one2)),
      path_cache(path_cache) {

    samebasis = false;

    calculate(config);
}

void HamiltonianTwo::calculate(const Configuration &conf_tot) {
    boost::filesystem::path path_cache_mat;
    if (utils::is_complex<scalar_t>::value) {
        path_cache_mat = path_cache / "cache_matrix_complex";
    } else {
        path_cache_mat = path_cache / "cache_matrix_real";
    }
    if (!boost::filesystem::exists(path_cache_mat)) {
        boost::filesystem::create_directory(path_cache_mat);
    }

    double tol = 1e-32;

    if (hamiltonian_one1->size() != hamiltonian_one2->size()) {
        throw std::runtime_error(
            "The number of single atom Hamiltonians must be the same for both atoms.");
    }

    size_t nSteps_one = hamiltonian_one1->size();

    // --- generate configuration ---
    std::vector<Configuration> conf_mat;
    conf_mat.reserve(nSteps_one);

    // new, pair hamiltonian specific configuration

    if (samebasis) {
        basis = std::make_shared<BasisnamesTwo>(hamiltonian_one1->names()); // TODO remove
    } else {
        basis = std::make_shared<BasisnamesTwo>(hamiltonian_one1->names(),
                                                hamiltonian_one2->names()); // TODO remove
    }
    Configuration conf_matpair = basis->getConf();
    conf_matpair["deltaEPair"] = conf_tot["deltaEPair"];
    conf_matpair["deltaNPair"] = conf_tot["deltaNPair"];
    conf_matpair["deltaLPair"] = conf_tot["deltaLPair"];
    conf_matpair["deltaJPair"] = conf_tot["deltaJPair"];
    conf_matpair["deltaMPair"] = conf_tot["deltaMPair"];
    // conf_matpair["conserveM"] = conf_tot["conserveM"];
    // conf_matpair["conserveParityL"] = conf_tot["conserveParityL"];
    conf_matpair["exponent"] = conf_tot["exponent"];

    for (size_t i = 0; i < nSteps_one; ++i) {
        // old, single atom hamiltonian specific configuration
        Configuration conf_matsingle = *hamiltonian_one1->getParams(i); // TODO
        conf_matsingle += conf_matpair;
        conf_mat.push_back(conf_matsingle);
        // conf_mat.push_back(conf_matsingle + conf_matpair); // conf_matpair overwrites settings in
        // conf_matsingle // TODO
    }

    // setup variables
    conf_mat.back()["species1"] >>
        species1; // TODO order state inside cinfiguration object config.order()
    conf_mat.back()["species2"] >> species2; // TODO order state inside cinfiguration object
    conf_mat.back()["deltaEPair"] >> deltaE;
    conf_mat.back()["deltaNPair"] >> deltaN;
    conf_mat.back()["deltaLPair"] >> deltaL;
    conf_mat.back()["deltaJPair"] >> deltaJ;
    conf_mat.back()["deltaMPair"] >> deltaM;
    // conserveM = conf_tot["conserveM"].str() == "true";
    // conserveParityL = conf_tot["conserveParityL"].str() == "true";
    conf_tot["steps"] >> nSteps_two;
    conf_tot["minR"] >> min_R;
    conf_tot["maxR"] >> max_R;
    conf_tot["exponent"] >> multipoleexponent;

    double minEx, minEy, minEz, maxEx, maxEy, maxEz, minBx, minBy, minBz, maxBx, maxBy, maxBz;
    conf_tot["minEx"] >> minEx;
    conf_tot["minEy"] >> minEy;
    conf_tot["minEz"] >> minEz;
    conf_tot["maxEx"] >> maxEx;
    conf_tot["maxEy"] >> maxEy;
    conf_tot["maxEz"] >> maxEz;
    conf_tot["minBx"] >> minBx;
    conf_tot["minBy"] >> minBy;
    conf_tot["minBz"] >> minBz;
    conf_tot["maxBx"] >> maxBx;
    conf_tot["maxBy"] >> maxBy;
    conf_tot["maxBz"] >> maxBz;
    // bool fields_change_m = (minEx != 0) || (minEy != 0) || (maxEx != 0) || (maxEy != 0) || (minBx
    // != 0) || (minBy != 0) || (maxBx != 0) || (maxBy != 0); // TODO wie richtig? so ist es eine
    // variable, die von mehreren matrizen abhaengt bool fields_change_l = (minEx != 0) || (minEy !=
    // 0) || (minEz != 0) || (maxEx != 0) || (maxEy != 0) || (maxEz != 0); // TODO wie richtig? so
    // ist es eine variable, die von mehreren matrizen abhaengt

    // fields_change_m = true; // TODO
    // fields_change_l = true; // TODO

    std::vector<parity_t> sym_inversion;
    if (conf_tot["invE"].str() == "true") {
        sym_inversion.push_back(EVEN);
    }
    if (conf_tot["invO"].str() == "true") {
        sym_inversion.push_back(ODD);
    }
    if (sym_inversion.empty()) {
        sym_inversion.push_back(NA);
    }

    std::vector<parity_t> sym_permutation;
    if (conf_tot["perE"].str() == "true") {
        sym_permutation.push_back(EVEN);
    }
    if (conf_tot["perO"].str() == "true") {
        sym_permutation.push_back(ODD);
    }
    if (sym_permutation.empty()) {
        sym_permutation.push_back(NA);
    }

    std::vector<parity_t> sym_reflection;
    if (conf_tot["refE"].str() == "true") {
        sym_reflection.push_back(EVEN);
    }
    if (conf_tot["refO"].str() == "true") {
        sym_reflection.push_back(ODD);
    }
    if (sym_reflection.empty()) {
        sym_reflection.push_back(NA);
    }

    bool conserveM = conf_tot["conserveM"].str() == "true";

    bool sametrafo = conf_tot["sametrafo"].str() == "true";

    bool zerotheta = conf_tot["zerotheta"].str() == "true";

    if (min_R == max_R && nSteps_one == 1) {
        nSteps_two = 1;
    }

    ////////////////////////////////////////////////////////
    ////// Restrict single atom states /////////////////////
    ////////////////////////////////////////////////////////

    // === Restrict states of atom 1 ===

    auto basis_one1 = hamiltonian_one1->names();
    std::vector<StateOneOld> initial1 = basis_one1->initial();
    std::vector<bool> necessary1(basis_one1->size(), false);

    for (const auto &state : *basis_one1) {
        bool validN = false;
        bool validL = false;
        bool validJ = false;
        bool validM = false;

        for (const auto &initial : initial1) {
            if (deltaN < 0 || std::abs(state.n - initial.n) <= deltaN) {
                validN = true;
            }
            if (deltaL < 0 || std::abs(state.l - initial.l) <= deltaL) {
                validL = true;
            }
            if (deltaJ < 0 || std::abs(state.j - initial.j) <= deltaJ) {
                validJ = true;
            }
            if (deltaM < 0 || std::abs(state.m - initial.m) <= deltaM) {
                validM = true;
            }
        }

        if (validN && validL && validJ && validM) {
            necessary1[state.idx] = true;
        }
    }

    hamiltonian_one1->removeUnnecessaryStates(necessary1);

    // === Restrict states of atom 2 ===

    if (!samebasis) {
        auto basis_one2 = hamiltonian_one2->names();
        std::vector<StateOneOld> initial2 = basis_one2->initial();
        std::vector<bool> necessary2(basis_one2->size(), false);

        for (const auto &state : *basis_one2) {
            bool validN = false;
            bool validL = false;
            bool validJ = false;
            bool validM = false;

            for (const auto &initial : initial2) {
                if (deltaN < 0 || std::abs(state.n - initial.n) <= deltaN) {
                    validN = true;
                }
                if (deltaL < 0 || std::abs(state.l - initial.l) <= deltaL) {
                    validL = true;
                }
                if (deltaJ < 0 || std::abs(state.j - initial.j) <= deltaJ) {
                    validJ = true;
                }
                if (deltaM < 0 || std::abs(state.m - initial.m) <= deltaM) {
                    validM = true;
                }
            }

            if (validN && validL && validJ && validM) {
                necessary2[state.idx] = true;
            }
        }

        hamiltonian_one2->removeUnnecessaryStates(necessary2);
    }

    ////////////////////////////////////////////////////////
    ////// Build pair state basis //////////////////////////
    ////////////////////////////////////////////////////////

    // === Build pair state basis ===

    std::cout << "Two-atom Hamiltonian, build pair state basis" << std::endl;

    if (samebasis) {
        basis = std::make_shared<BasisnamesTwo>(hamiltonian_one1->names());
    } else {
        basis =
            std::make_shared<BasisnamesTwo>(hamiltonian_one1->names(), hamiltonian_one2->names());
    }

    std::cout << "Two-atom Hamiltonian, basis size without restrictions: " << basis->size()
              << std::endl;

    // === Determine necessary symmetries ===
    std::cout << "Two-atom Hamiltonian, determine symmetrized subspaces" << std::endl;

    StateTwoOld initial = basis->initial();
    parity_t initalParityL = (std::pow(-1, initial.l[0] + initial.l[1]) > 0) ? EVEN : ODD;
    int initalM = initial.m[0] + initial.m[1];
    int initalJ = initial.j[0] + initial.j[1];
    bool samestates = initial.first() == initial.second();

    std::vector<int> sym_rotation;
    if (conserveM) {
        if (!sametrafo) {
            for (const auto &state : *basis) {
                sym_rotation.push_back(state.m[0] + state.m[1]);
            }
        } else if (!zerotheta) {
            for (int M = -initalJ; M <= initalJ; ++M) {
                sym_rotation.push_back(M);
            }
        } else {
            sym_rotation.push_back(initalM);
        }
    } else {
        sym_rotation.push_back(NA);
    }

    Symmetry sym;
    std::set<Symmetry> symmetries_set;
    for (const parity_t &inv : sym_inversion) {
        sym.inversion = inv;

        // In case of even inversion symmetry and the same inital state for the first and second
        // atom: the inital state ist not contained in the corresponding block
        if (sym.inversion == EVEN && samestates && sametrafo) {
            continue;
        }

        for (const parity_t &per : sym_permutation) {
            sym.permutation = per;

            // In case of even permutation symmetry and the same inital state for the first and
            // second atom: the inital state ist not contained in the corresponding block
            if (sym.permutation == EVEN && samestates && sametrafo) {
                continue;
            }

            // In case of inversion and permutation symmetry: the orbital parity is conserved and we
            // just use the blocks with the same parity as the inital state
            if (sym.inversion != NA && sym.permutation != NA && sametrafo) {
                if (sym.inversion * sym.permutation != initalParityL) {
                    continue;
                }
            }

            for (const parity_t &ref : sym_reflection) {
                sym.reflection = ref;
                for (const int &rot : sym_rotation) {

                    // In case of reflection symmetry: do just use the absolute value of rot
                    if (sym.reflection != NA) {
                        sym.rotation = std::abs(rot);
                    } else {
                        sym.rotation = rot;
                    }

                    symmetries_set.insert(sym);
                }
            }
        }
    }
    std::vector<Symmetry> symmetries(symmetries_set.begin(), symmetries_set.end());

    // === Build up the list of necessary pair states ===
    std::cout << "Two-atom Hamiltonian, build up the list of necessary pair states" << std::endl;

    // Apply energy cutoff
    std::vector<bool> necessary_tmp(basis->size(), false);

    auto nSteps_one_i = static_cast<int>(nSteps_one);

#pragma omp parallel for
    for (int i = 0; i < nSteps_one_i; ++i) {
        energycutoff(*(hamiltonian_one1->get(i)), *(hamiltonian_one2->get(i)), deltaE,
                     necessary_tmp);
    }

    // Apply restrictions due to symmetries
    std::vector<bool> necessary(basis->size(), false);

    for (const auto &state : *basis) {
        for (Symmetry sym : symmetries) {
            float M = state.m[0] + state.m[1];
            int parityL = std::pow(-1, state.l[0] + state.l[1]);

            // In case of rotation symmetry: skip pair states with wrong total magnetic momentum
            if (sym.rotation != NA && sym.rotation != M &&
                (sym.reflection == NA || sym.rotation != -M)) {
                continue;
            }

            // In case of inversion and permutation symmetry: skip pair states with wrong orbital
            // parity
            if (sym.inversion != NA && sym.permutation != NA && parityL != initalParityL &&
                sametrafo) {
                continue;
            }

            necessary[state.idx] = necessary_tmp[state.idx];
        }
    }

    int numNecessary = std::count(necessary.begin(), necessary.end(), true);
    std::cout << "Two-atom Hamiltonian, basis size with restrictions: " << numNecessary
              << std::endl;
    std::cout << boost::format(">>BAS%7d") % numNecessary << std::endl;

    // === Save pair state basis ===
    std::cout << "Two-atom Hamiltonian, save pair state basis" << std::endl;

    // initialize uuid generator
    boost::uuids::random_generator generator;

    // generate uuid
    std::string uuid;
    boost::uuids::uuid u = generator();
    boost::algorithm::hex(u.begin(), u.end(), std::back_inserter(uuid));

    // save pair state basis
    boost::filesystem::path path_basis = boost::filesystem::temp_directory_path();
    path_basis /= "basis_two_" + uuid + ".csv";
    basis->save(
        path_basis
            .string()); // TODO save only necessary entries, i.e. save pair state basis in sparse
                        // format (possibility, remove basis states but keep their idx - this would
                        // also make "if (necessary) continue" unneeded; then, "combine" has to
                        // check existence of basis element and the python script has to be adapted)

    std::cout << boost::format(">>STA %s") % path_basis.string() << std::endl;

    ////////////////////////////////////////////////////////
    ////// Construct atom-atom interaction /////////////////
    ////////////////////////////////////////////////////////

    // Construct pair Hamiltonians for all orders of the multipole expansion

    std::vector<int> exponent_multipole;
    std::vector<Hamiltonianmatrix> mat_multipole;
    MatrixElements matrixelements_atom1(conf_tot, species1,
                                        (path_cache / "cache_elements.db").string());
    MatrixElements matrixelements_atom2(conf_tot, species2,
                                        (path_cache / "cache_elements.db").string());
    std::vector<idx_t> size_mat_multipole;

    int idx_multipole_max = -1;

    if (multipoleexponent > 2) {

        // --- Initialize two-atom interaction Hamiltonians ---
        std::cout << "Two-atom Hamiltonian, initialize interaction Hamiltonians" << std::endl;

        int kappa_min = 1; // spherical dipole operators
        int kappa_max = multipoleexponent - kappa_min - 1;
        int sumOfKappas_min = kappa_min + kappa_min;
        int sumOfKappas_max = kappa_max + kappa_min;
        idx_multipole_max = sumOfKappas_max - sumOfKappas_min;

        exponent_multipole.reserve(idx_multipole_max + 1);
        mat_multipole.reserve(idx_multipole_max + 1);
        size_mat_multipole.resize(idx_multipole_max + 1);

        // --- Precalculate matrix elements --- // TODO parallelization
        std::cout << "Two-atom Hamiltonian, get one-atom states needed for the pair state basis"
                  << std::endl;

        auto basis_one1_needed = std::make_shared<BasisnamesOne>(BasisnamesOne::fromFirst(basis));
        auto basis_one2_needed = std::make_shared<BasisnamesOne>(BasisnamesOne::fromSecond(basis));

        for (int kappa = kappa_min; kappa <= kappa_max; ++kappa) {
            std::cout << "Two-atom Hamiltonian, precalculate matrix elements for kappa = " << kappa
                      << std::endl;
            matrixelements_atom1.precalculateMultipole(basis_one1_needed, kappa);
            matrixelements_atom2.precalculateMultipole(basis_one2_needed, kappa);
        }

        // TODO if (samebasis) ...

        // --- Count entries of two-atom interaction Hamiltonians ---
        std::cout
            << "Two-atom Hamiltonian, count number of entries within the interaction Hamiltonians"
            << std::endl;

        for (int sumOfKappas = sumOfKappas_min; sumOfKappas <= sumOfKappas_max; ++sumOfKappas) {
            int idx_multipole = sumOfKappas - sumOfKappas_min;

            for (const auto &state_col : *basis) { // TODO parallelization
                if (!necessary[state_col.idx]) {
                    continue;
                }

                int M_col = state_col.first().m + state_col.second().m;

                for (const auto &state_row : *basis) {
                    if (!necessary[state_row.idx]) {
                        continue;
                    }

                    if (state_row.idx < state_col.idx) {
                        continue;
                    }
                    int M_row = state_row.first().m + state_row.second().m;
                    if (M_col != M_row) {
                        continue;
                    }

                    // multipole interaction with 1/R^(sumOfKappas+1) = 1/R^(idx_multipole+3) decay
                    for (int kappa1 = kappa_min; kappa1 <= sumOfKappas - 1; ++kappa1) {
                        int kappa2 = sumOfKappas - kappa1;

                        // allowed deltaL, deltaJ, and deltaM?
                        if (selectionRulesMultipole(state_row.first(), state_col.first(), kappa1) &&
                            selectionRulesMultipole(state_row.second(), state_col.second(),
                                                    kappa2)) {
                            int q1 = state_row.first().m - state_col.first().m;
                            int q2 = state_row.second().m - state_col.second().m;

                            // total momentum preserved?
                            if (q1 == -q2) {
                                size_mat_multipole[idx_multipole]++;
                                break;
                            }
                        }
                    }
                }
            }
        }

        // --- Construct two-atom interaction Hamiltonians ---
        size_t size_basis = basis->size();

        for (int sumOfKappas = sumOfKappas_min; sumOfKappas <= sumOfKappas_max; ++sumOfKappas) {
            std::cout
                << "Two-atom Hamiltonian, construct interaction Hamiltonian that belongs to 1/R^"
                << sumOfKappas + 1 << std::endl;

            int idx_multipole = sumOfKappas - sumOfKappas_min;

            exponent_multipole.push_back(sumOfKappas + 1);
            mat_multipole.emplace_back(
                size_basis,
                2 * size_mat_multipole[idx_multipole]); // factor of 2 because triangular matrix is
                                                        // not sufficient

            for (const auto &state_col : *basis) { // TODO parallelization
                if (!necessary[state_col.idx]) {
                    continue;
                }

                int M_col = state_col.first().m + state_col.second().m;

                for (const auto &state_row : *basis) {
                    if (!necessary[state_row.idx]) {
                        continue;
                    }

                    if (state_row.idx < state_col.idx) {
                        continue;
                    }
                    int M_row = state_row.first().m + state_row.second().m;
                    if (M_col != M_row) {
                        continue;
                    }

                    // construct basis
                    if (state_row.idx == state_col.idx) {
                        mat_multipole[idx_multipole].addBasis(state_row.idx, state_col.idx, 1);
                    }

                    // multipole interaction with 1/R^(sumOfKappas+1) = 1/R^(idx_multipole+3) decay
                    double val = 0;

                    for (int kappa1 = kappa_min; kappa1 <= sumOfKappas - 1; ++kappa1) {
                        int kappa2 = sumOfKappas - kappa1;

                        // allowed deltaL, deltaJ, and deltaM?
                        if (selectionRulesMultipole(state_row.first(), state_col.first(), kappa1) &&
                            selectionRulesMultipole(state_row.second(), state_col.second(),
                                                    kappa2)) {
                            int q1 = state_row.first().m - state_col.first().m;
                            int q2 = state_row.second().m - state_col.second().m;

                            // total momentum preserved?
                            if (q1 == -q2) {
                                double binomials = boost::math::binomial_coefficient<double>(
                                                       kappa1 + kappa2, kappa1 + q1) *
                                    boost::math::binomial_coefficient<double>(kappa1 + kappa2,
                                                                              kappa2 - q2);
                                val += inverse_electric_constant * std::pow(-1, kappa2) *
                                    std::sqrt(binomials) *
                                    matrixelements_atom1.getMultipole(state_row.first(),
                                                                      state_col.first(), kappa1) *
                                    matrixelements_atom2.getMultipole(state_row.second(),
                                                                      state_col.second(), kappa2);
                            }
                        }
                    }

                    if (std::abs(val) > tol) {
                        mat_multipole[idx_multipole].addEntries(state_row.idx, state_col.idx, val);
                        if (state_row.idx != state_col.idx) {
                            mat_multipole[idx_multipole].addEntries(
                                state_col.idx, state_row.idx,
                                val); // triangular matrix is not sufficient because of basis change
                        }
                    }

                    // TODO state_two soll std::array<state_one, 2> sein! Dann geht auch die Abfrage
                    // der selection rules eindeutiger
                }
            }

            std::cout
                << "Two-atom Hamiltonian, compress interaction Hamiltonian that belongs to 1/R^"
                << sumOfKappas + 1 << std::endl;

            mat_multipole[idx_multipole].compress(basis->dim(),
                                                  basis->dim()); // TODO substitute dim() by size()
        }
    }

    ////////////////////////////////////////////////////////
    ////// Prepare processing of Hamiltonians //////////////
    ////////////////////////////////////////////////////////

    // TODO Put the logic in its own class

    std::cout << "Two-atom Hamiltonian, process Hamiltonians" << std::endl;

    // === Open database ===
    boost::filesystem::path path_db;

    if (utils::is_complex<scalar_t>::value) {
        path_db = path_cache / "cache_matrix_complex.db";
    } else {
        path_db = path_cache / "cache_matrix_real.db";
    }
    sqlite::handle db(path_db.string());
    sqlite::statement stmt(db);

    // === Initialize variables ===
    bool flag_perhapsmissingtable = true;

    matrix_path.resize(nSteps_two * symmetries.size());

    auto indices_symmetry_i = static_cast<int>(symmetries.size());

    // --- Determine combined single atom matrices ---
    // Construct pair Hamiltonian consistent of combined one-atom Hamiltonians (1 x Hamiltonian2 +
    // Hamiltonian1 x 1)

    std::vector<Hamiltonianmatrix> mat_single;

    // Check if one_atom Hamiltonians change with step_two
    // It is assumed that nSteps_one = 1 if nSteps_two != nSteps_one // TODO introduce variable
    // "is_mat_single_const" to improve readability
    if (nSteps_two != nSteps_one) {
        std::cout
            << "Two-atom Hamiltonian, construct contribution of combined one-atom Hamiltonians"
            << std::endl;

        mat_single.resize(symmetries.size());

#pragma omp parallel for
        for (int idx_symmetry = 0; idx_symmetry < indices_symmetry_i; ++idx_symmetry) {
            Symmetry sym = symmetries[idx_symmetry];

            // Combine the Hamiltonians of the two atoms
            mat_single[idx_symmetry] = combine(*(hamiltonian_one1->get(0)),
                                               *(hamiltonian_one2->get(0)), deltaE, basis, sym);

            // Remove more or less empty basis vectors
            mat_single[idx_symmetry].removeUnnecessaryBasisvectors();
        }
    }

    // --- Determine transformed interaction matrices ---
    std::vector<Hamiltonianmatrix> mat_multipole_transformed;

    // Check if one_atom Hamiltonians change with step_two
    if (nSteps_two != nSteps_one) {
        std::cout << "Two-atom Hamiltonian, construct transformed interaction matrices"
                  << std::endl;

        mat_multipole_transformed.resize(symmetries.size() * (idx_multipole_max + 1));

#pragma omp parallel for
        for (int idx_symmetry = 0; idx_symmetry < indices_symmetry_i; ++idx_symmetry) {
            for (int idx_multipole = 0; idx_multipole <= idx_multipole_max; ++idx_multipole) {
                mat_multipole_transformed[idx_symmetry * (idx_multipole_max + 1) + idx_multipole] =
                    mat_multipole[idx_multipole].changeBasis(mat_single[idx_symmetry].basis());
            }
        }
    }

    ////////////////////////////////////////////////////////
    ////// Loop through steps and symmetries ///////////////
    ////////////////////////////////////////////////////////

    std::cout << boost::format(">>TOT%7d") % (nSteps_two * symmetries.size()) << std::endl;

    auto nSteps_two_i = static_cast<int>(nSteps_two);

#pragma omp parallel for schedule(static, 1)

    // Loop through steps
    for (int step_two = 0; step_two < nSteps_two_i; ++step_two) {

        // Loop through symmetries
        for (size_t idx_symmetry = 0; idx_symmetry < symmetries.size(); ++idx_symmetry) {
            Symmetry sym = symmetries[idx_symmetry];

            size_t step = step_two * symmetries.size() + idx_symmetry;

            // === Get parameters for the current position inside the loop ===
            int single_idx = (nSteps_two == nSteps_one) ? step_two : 0;

            // Get interatomic distance
            double normalized_position = (nSteps_two > 1) ? step_two / (nSteps_two - 1.) : 0;
            double position = min_R + normalized_position * (max_R - min_R);

            // Get configuration and save postions and symmetries
            Configuration conf = conf_mat[single_idx];
            conf["R"] << position;
            conf["inversion"] << ((sym.inversion == NA) ? "" : std::to_string(sym.inversion));
            conf["permutation"] << ((sym.permutation == NA) ? "" : std::to_string(sym.permutation));
            conf["reflection"] << ((sym.reflection == NA) ? "" : std::to_string(sym.reflection));
            conf["rotation"] << ((sym.rotation == NA) ? "" : std::to_string(sym.rotation));

            // === Create table if necessary ===
            std::stringstream query;
            std::string spacer;

            if (flag_perhapsmissingtable) {
                query << "CREATE TABLE IF NOT EXISTS cache_two (uuid text NOT NULL PRIMARY KEY, "
                         "created TIMESTAMP DEFAULT CURRENT_TIMESTAMP, "
                         "accessed TIMESTAMP DEFAULT CURRENT_TIMESTAMP";
                for (auto p : conf) {
                    query << ", " << p.first << " text";
                }
                query << ", UNIQUE (";
                for (auto p : conf) {
                    query << spacer << p.first;
                    spacer = ", ";
                }
                query << "));";

#pragma omp critical(database)
                stmt.exec(query.str());

                flag_perhapsmissingtable = false;
            }

            // === Get uuid as filename ===
            std::string uuid;
            query.str(std::string());
            spacer = "";
            query << "SELECT uuid FROM cache_two WHERE ";
            for (auto p : conf) {
                query << spacer << p.first << "='" << p.second.str() << "'";
                spacer = " AND ";
            }
            query << ";";

#pragma omp critical(database)
            {
                sqlite::statement stmt(db, query.str());
                stmt.prepare();
                if (stmt.step()) {
                    uuid = stmt.get<std::string>(0);
                }
            }

            if (!uuid.empty()) {
                query.str(std::string());
                query << "UPDATE cache_two SET accessed = CURRENT_TIMESTAMP WHERE uuid = '" << uuid
                      << "';";
#pragma omp critical(database)
                stmt.exec(query.str()); // TODO check whether this slows down the program

            } else {
                boost::uuids::uuid u = generator();
                boost::algorithm::hex(u.begin(), u.end(), std::back_inserter(uuid));

                query.str(std::string());
                query << "INSERT INTO cache_two (uuid";
                for (auto p : conf) {
                    query << ", " << p.first;
                }
                query << ") values ( '" << uuid << "'";
                for (auto p : conf) {
                    query << ", "
                          << "'" << p.second.str() << "'";
                }
                query << ");";
#pragma omp critical(database)
                stmt.exec(query.str());
            }

            // === Check existence of files ===

            // Check whether .mat and .json file exists and compare settings in program with
            // settings in .json file
            boost::filesystem::path path, path_mat, path_json;

            path = path_cache_mat / ("two_" + uuid);
            path_mat = path;
            path_mat.replace_extension(".mat");
            path_json = path;
            path_json.replace_extension(".json");

            bool is_existing = false;
            if (boost::filesystem::exists(path_mat)) {
                if (boost::filesystem::exists(path_json)) {
                    Configuration params_loaded;
                    params_loaded.load_from_json(path_json.string());
                    if (conf == params_loaded) {
                        is_existing = true;
                    }
                }
            }

            // Create .json file if "is_existing" is false
            if (!is_existing) {
                conf.save_to_json(path_json.string());
            }

            // === Build and diagonalize total matrix if not existent ===
            Hamiltonianmatrix totalmatrix;

            if (!is_existing || !totalmatrix.load(path_mat.string())) {

                // --- Combine single atom matrices ---
                if (nSteps_two == nSteps_one) {
                    totalmatrix = combine(*(hamiltonian_one1->get(step_two)),
                                          *(hamiltonian_one2->get(step_two)), deltaE, basis, sym);
                    totalmatrix.removeUnnecessaryBasisvectors();
                } else {
                    totalmatrix = mat_single[idx_symmetry];
                }

                // --- Add interaction ---
                for (int idx_multipole = 0; idx_multipole <= idx_multipole_max; ++idx_multipole) {
                    double pos = 1. / std::pow(position, exponent_multipole[idx_multipole]);
                    if (nSteps_two == nSteps_one) {
                        totalmatrix +=
                            mat_multipole[idx_multipole].changeBasis(totalmatrix.basis()) * pos;
                    } else {
                        totalmatrix +=
                            mat_multipole_transformed[idx_symmetry * (idx_multipole_max + 1) +
                                                      idx_multipole] *
                            pos;
                    }
                }

                // Stdout: Hamiltonian assembled
#pragma omp critical(textoutput)
                {
                    std::cout << boost::format(">>DIM%7d") % totalmatrix.num_basisvectors()
                              << std::endl;
                    std::cout << "Two-atom Hamiltonian, " << step + 1 << ". Hamiltonian assembled"
                              << std::endl;
                }

                // --- Diagonalize matrix and save diagonalized matrix ---
                totalmatrix.diagonalize();
                totalmatrix.save(path_mat.string());

                // Stdout: Hamiltonian diagonalized
#pragma omp critical(textoutput)
                {
                    std::cout << boost::format(">>OUT%7d%7d%7d%7d %s") % (step + 1) % step_two %
                            symmetries.size() % idx_symmetry % path.string()
                              << std::endl;
                    std::cout << "Two-atom Hamiltonian, " << step + 1
                              << ". Hamiltonian diagonalized" << std::endl;
                }
            } else {

                // Stdout: Hamiltonian loaded
#pragma omp critical(textoutput)
                {
                    std::cout << boost::format(">>DIM%7d") % totalmatrix.num_basisvectors()
                              << std::endl;
                    std::cout << boost::format(">>OUT%7d%7d%7d%7d %s") % (step + 1) % step_two %
                            symmetries.size() % idx_symmetry % path.string()
                              << std::endl;
                    std::cout << "Two-atom Hamiltonian, " << step + 1 << ". Hamiltonian loaded"
                              << std::endl;
                }
            }

            // === Store path to configuration and diagonalized matrix ===
            matrix_path[step] = path.string();
        }
    }

    std::cout << "Two-atom Hamiltonian, all Hamiltonians processed" << std::endl;
}
