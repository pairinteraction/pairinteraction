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

#include "communication.h"
#include "HamiltonianOne.h"
#include <stdexcept>

HamiltonianOne::HamiltonianOne(const Configuration &config, boost::filesystem::path& path_cache, std::shared_ptr<BasisnamesOne> basis_one) : Hamiltonian<BasisnamesOne>(), path_cache(path_cache) {
    basis = basis_one;
    configure(config);
    build();
}

const Configuration& HamiltonianOne::getConf() const { // TODO in Configurable Klasse auslagern, von der geerbt werrden soll
    return basicconf;
}

void HamiltonianOne::changeToSpherical(real_t val_x, real_t val_y, real_t val_z, real_t& val_p, real_t& val_m, real_t& val_0) {
    if(val_y != 0) {
        std::string msg("For fields with non-zero y-coordinates, a complex data type is needed.");
        auto context = zmq::context();
        auto publisher = context.socket(ZMQ_PUB);
        publisher.connect("tcp://localhost:5556");
        publisher.send(">>ERR%s", msg.c_str());
        throw std::runtime_error(msg);
    }
    val_p = -val_x/std::sqrt(2);
    val_m = val_x/std::sqrt(2);
    val_0 = val_z;
}

void HamiltonianOne::changeToSpherical(real_t val_x, real_t val_y, real_t val_z, std::complex<real_t>& val_p, std::complex<real_t>& val_m, std::complex<real_t>& val_0) {
    val_p = std::complex<real_t>(-val_x/std::sqrt(2),-val_y/std::sqrt(2));
    val_m = std::complex<real_t>(val_x/std::sqrt(2),-val_y/std::sqrt(2));
    val_0 = std::complex<real_t>(val_z,0);
}

void HamiltonianOne::configure(const Configuration &config) {
    basicconf = basis->getConf();
    basicconf["deltaESingle"] << config["deltaESingle"];
    basicconf["diamagnetism"] << config["diamagnetism"];

    basicconf["deltaESingle"] >> deltaE;
    basicconf["species1"] >> species;

    diamagnetism = basicconf["diamagnetism"].str() == "true";

    config["minBx"] >> min_B_x;
    config["minBy"] >> min_B_y;
    config["minBz"] >> min_B_z;
    config["minEx"] >> min_E_x;
    config["minEy"] >> min_E_y;
    config["minEz"] >> min_E_z;
    config["maxBx"] >> max_B_x;
    config["maxBy"] >> max_B_y;
    config["maxBz"] >> max_B_z;
    config["maxEx"] >> max_E_x;
    config["maxEy"] >> max_E_y;
    config["maxEz"] >> max_E_z;

    if ((min_B_x == max_B_x) &&
            (min_B_y == max_B_y) &&
            (min_B_z == max_B_z) &&
            (min_E_x == max_E_x) &&
            (min_E_y == max_E_y) &&
            (min_E_z == max_E_z)) {
        nSteps = 1;
    } else {
        config["steps"] >> nSteps;
    }
}

void HamiltonianOne::build() {
    boost::filesystem::path path_cache_mat;
    if (utils::is_complex<scalar_t>::value) {
        path_cache_mat = path_cache / "cache_matrix_complex";
    } else {
        path_cache_mat = path_cache / "cache_matrix_real";
    }
    if(!boost::filesystem::exists(path_cache_mat)){
        boost::filesystem::create_directory(path_cache_mat);
    }

    real_t tol = 1e-32;

    ////////////////////////////////////////////////////////
    ////// Build single atom basis and Hamiltonian /////////
    ////////////////////////////////////////////////////////

    // === Calculate one-atom Hamiltonian ===

    // --- Count entries of one-atom Hamiltonian ---
    size_t size_basis = basis->size();
    size_t size_energy = basis->size();

    // --- Construct one-atom  Hamiltonian and basis ---
    std::cout << "One-atom Hamiltonian, construct diagonal Hamiltonian" << std::endl;

    Hamiltonianmatrix hamiltonian_energy(size_basis, size_energy);

    real_t energy_initial = 0;
    for (const auto &state: basis->initial()) {
        energy_initial += energy_level(species,state.n,state.l,state.j);
    }
    energy_initial /= basis->initial().size(); // TODO save it to the json file

    std::vector<bool> is_necessary(basis->size(),false);
    idx_t idx = 0;
    for (const auto &state : *basis) {
        real_t val = energy_level(species,state.n,state.l,state.j)-energy_initial;
        if (std::abs(val) < deltaE+1e-11 || deltaE < 0) { // TODO
            is_necessary[state.idx] = true;
            hamiltonian_energy.addEntries(idx,idx,val);
            hamiltonian_energy.addBasis(idx,idx,1);
            ++idx;
        }
    }
    std::cout << "One-atom Hamiltonian, basis size without restrictions: " << basis->size() << std::endl;

    basis->removeUnnecessaryStates(is_necessary);

    hamiltonian_energy.compress(basis->dim(), basis->dim());

    auto context = zmq::context();
    auto publisher = context.socket(ZMQ_PUB);
    publisher.connect("tcp://localhost:5556");

    std::cout << "One-atom Hamiltonian, basis size with restrictions: " << basis->size() << std::endl;
    publisher.send(">>BAS%7d", basis->size());

    // === Save single atom basis ===
    std::cout << "One-atom Hamiltonian, save single atom basis" << std::endl;

    // initialize uuid generator
    boost::uuids::random_generator generator;

    // generate uuid
    std::string uuid;
    boost::uuids::uuid u = generator();
    boost::algorithm::hex(u.begin(), u.end(), std::back_inserter(uuid));

    // save basis
    boost::filesystem::path path_basis = boost::filesystem::temp_directory_path();
    path_basis /= "basis_one_"+uuid+".csv";
    basis->save(path_basis.string());

    publisher.send(">>STA %s", path_basis.c_str());


    ////////////////////////////////////////////////////////
    ////// Construct atom-field interaction ////////////////
    ////////////////////////////////////////////////////////

    // --- Obtain existence of fields ---
    scalar_t min_E_0, min_E_p, min_E_m, min_B_0, min_B_p, min_B_m, max_E_0, max_E_p, max_E_m, max_B_0, max_B_p, max_B_m;
    changeToSpherical(min_E_x, min_E_y, min_E_z, min_E_p, min_E_m, min_E_0);
    changeToSpherical(max_E_x, max_E_y, max_E_z, max_E_p, max_E_m, max_E_0);
    changeToSpherical(min_B_x, min_B_y, min_B_z, min_B_p, min_B_m, min_B_0);
    changeToSpherical(max_B_x, max_B_y, max_B_z, max_B_p, max_B_m, max_B_0);

    bool exist_E_0 = (std::abs(min_E_0) != 0 || std::abs(max_E_0) != 0);
    bool exist_E_1 = (std::abs(min_E_p) != 0 || std::abs(max_E_p) != 0);
    bool exist_B_0 = (std::abs(min_B_0) != 0 || std::abs(max_B_0) != 0);
    bool exist_B_1 = (std::abs(min_B_p) != 0 || std::abs(max_B_p) != 0);

    // --- Precalculate matrix elements --- // TODO parallelization
    std::cout << "One-atom Hamiltonian, precalculate matrix elements" << std::endl;

    MatrixElements matrix_elements(basicconf, species, (path_cache / "cache_elements.db").string());

    if (exist_E_0) matrix_elements.precalculateElectricMomentum(basis, 0);
    if (exist_E_1) matrix_elements.precalculateElectricMomentum(basis, 1);
    if (exist_E_1) matrix_elements.precalculateElectricMomentum(basis, -1);

    if (exist_B_0) matrix_elements.precalculateMagneticMomentum(basis, 0);
    if (exist_B_1) matrix_elements.precalculateMagneticMomentum(basis, 1);
    if (exist_B_1) matrix_elements.precalculateMagneticMomentum(basis, -1);

    if (diamagnetism && (exist_B_0 || exist_B_1)) matrix_elements.precalculateDiamagnetism(basis, 0, 0);
    if (diamagnetism && (exist_B_0 || exist_B_1)) matrix_elements.precalculateDiamagnetism(basis, 2, 0);
    if (diamagnetism && exist_B_0 && exist_B_1) matrix_elements.precalculateDiamagnetism(basis, 2, 1);
    if (diamagnetism && exist_B_0 && exist_B_1) matrix_elements.precalculateDiamagnetism(basis, 2, -1);
    if (diamagnetism && exist_B_1) matrix_elements.precalculateDiamagnetism(basis, 2, 2);
    if (diamagnetism && exist_B_1) matrix_elements.precalculateDiamagnetism(basis, 2, -2);

    // --- Count entries of atom-field Hamiltonian ---
    std::cout << "One-atom Hamiltonian, count number of entries within the field Hamiltonian" << std::endl;

    size_basis = basis->size();
    size_t size_electricMomentum_0 = 0;
    size_t size_electricMomentum_p = 0;
    size_t size_electricMomentum_m = 0;

    size_t size_magneticMomentum_0 = 0;
    size_t size_magneticMomentum_p = 0;
    size_t size_magneticMomentum_m = 0;

    size_t size_diamagnetism_00 = 0;
    size_t size_diamagnetism_20 = 0;
    size_t size_diamagnetism_2p = 0;
    size_t size_diamagnetism_2m = 0;
    size_t size_diamagnetism_2pp = 0;
    size_t size_diamagnetism_2mm = 0;

    for (const auto &state_col : *basis) { // TODO parallelization
        for (const auto &state_row : *basis) {
            if (state_row.idx < state_col.idx) { // lower triangle only
                continue;
            }

            if (exist_E_0 && selectionRulesMultipole(state_row, state_col, 1, 0) ) {
                size_electricMomentum_0++;
            } else if (exist_E_1 && selectionRulesMultipole(state_row, state_col, 1, 1) ) {
                size_electricMomentum_p++;
            } else if (exist_E_1 && selectionRulesMultipole(state_row, state_col, 1, -1) ) {
                size_electricMomentum_m++;
            }

            if (exist_B_0 &&  selectionRulesMomentum(state_row, state_col, 0) ) {
                size_magneticMomentum_0++;
            } else if (exist_B_1 && selectionRulesMomentum(state_row, state_col, 1) ) {
                size_magneticMomentum_p++;
            } else if (exist_B_1 && selectionRulesMomentum(state_row, state_col, -1) ) {
                size_magneticMomentum_m++;
            }

            if (diamagnetism && (exist_B_0 || exist_B_1) && selectionRulesMultipole(state_row, state_col, 0, 0)) {
                size_diamagnetism_00++;
            } else if (diamagnetism && (exist_B_0 || exist_B_1) && selectionRulesMultipole(state_row, state_col, 2, 0)) {
                size_diamagnetism_20++;
            } else if (diamagnetism && (exist_B_0 && exist_B_1) && selectionRulesMultipole(state_row, state_col, 2, 1)) {
                size_diamagnetism_2p++;
            } else if (diamagnetism && (exist_B_0 && exist_B_1) && selectionRulesMultipole(state_row, state_col, 2, -1)) {
                size_diamagnetism_2m++;
            } else if (diamagnetism && (exist_B_1) && selectionRulesMultipole(state_row, state_col, 2, 2)) {
                size_diamagnetism_2pp++;
            } else if (diamagnetism && (exist_B_1) && selectionRulesMultipole(state_row, state_col, 2, -2)) {
                size_diamagnetism_2mm++;
            }
        }
    }

    // --- Construct atom-field Hamiltonian ---
    std::cout << "One-atom Hamiltonian, construct field Hamiltonian" << std::endl;

    Hamiltonianmatrix hamiltonian_electricMomentum_0(size_basis, size_electricMomentum_0);
    Hamiltonianmatrix hamiltonian_electricMomentum_p(size_basis, size_electricMomentum_p);
    Hamiltonianmatrix hamiltonian_electricMomentum_m(size_basis, size_electricMomentum_m);

    Hamiltonianmatrix hamiltonian_magneticMomentum_0(size_basis, size_magneticMomentum_0);
    Hamiltonianmatrix hamiltonian_magneticMomentum_p(size_basis, size_magneticMomentum_p);
    Hamiltonianmatrix hamiltonian_magneticMomentum_m(size_basis, size_magneticMomentum_m);

    Hamiltonianmatrix hamiltonian_diamagnetism_00(size_basis, size_diamagnetism_00);
    Hamiltonianmatrix hamiltonian_diamagnetism_20(size_basis, size_diamagnetism_20);
    Hamiltonianmatrix hamiltonian_diamagnetism_2p(size_basis, size_diamagnetism_2p);
    Hamiltonianmatrix hamiltonian_diamagnetism_2m(size_basis, size_diamagnetism_2m);
    Hamiltonianmatrix hamiltonian_diamagnetism_2pp(size_basis, size_diamagnetism_2pp);
    Hamiltonianmatrix hamiltonian_diamagnetism_2mm(size_basis, size_diamagnetism_2mm);

    for (const auto &state_col : *basis) { // TODO parallelization
        for (const auto &state_row : *basis) {
            if (state_row.idx < state_col.idx) {
                continue;
            }

            if (state_row.idx == state_col.idx) {
                hamiltonian_electricMomentum_0.addBasis(state_row.idx,state_col.idx,1);
                hamiltonian_electricMomentum_p.addBasis(state_row.idx,state_col.idx,1);
                hamiltonian_electricMomentum_m.addBasis(state_row.idx,state_col.idx,1);

                hamiltonian_magneticMomentum_0.addBasis(state_row.idx,state_col.idx,1);
                hamiltonian_magneticMomentum_p.addBasis(state_row.idx,state_col.idx,1);
                hamiltonian_magneticMomentum_m.addBasis(state_row.idx,state_col.idx,1);

                hamiltonian_diamagnetism_00.addBasis(state_row.idx,state_col.idx,1);
                hamiltonian_diamagnetism_20.addBasis(state_row.idx,state_col.idx,1);
                hamiltonian_diamagnetism_2p.addBasis(state_row.idx,state_col.idx,1);
                hamiltonian_diamagnetism_2m.addBasis(state_row.idx,state_col.idx,1);
                hamiltonian_diamagnetism_2pp.addBasis(state_row.idx,state_col.idx,1);
                hamiltonian_diamagnetism_2mm.addBasis(state_row.idx,state_col.idx,1);
            }

            if (exist_E_0 && selectionRulesMultipole(state_row, state_col, 1, 0) ) {
                real_t val = matrix_elements.getElectricMomentum(state_row, state_col);
                if (std::abs(val) > tol) {
                    hamiltonian_electricMomentum_0.addEntries(state_row.idx,state_col.idx,val);
                }
            } else if (exist_E_1 && selectionRulesMultipole(state_row, state_col, 1, 1) ) {
                real_t val = matrix_elements.getElectricMomentum(state_row, state_col);
                if (std::abs(val) > tol) {
                    hamiltonian_electricMomentum_p.addEntries(state_row.idx,state_col.idx,val);
                }
            } else if (exist_E_1 && selectionRulesMultipole(state_row, state_col, 1, -1) ) {
                real_t val = matrix_elements.getElectricMomentum(state_row, state_col);
                if (std::abs(val) > tol) {
                    hamiltonian_electricMomentum_m.addEntries(state_row.idx,state_col.idx,val);
                }
            }

            if (exist_B_0 && selectionRulesMomentum(state_row, state_col, 0) ) {
                real_t val = matrix_elements.getMagneticMomentum(state_row, state_col);
                if (std::abs(val) > tol) {
                    hamiltonian_magneticMomentum_0.addEntries(state_row.idx,state_col.idx,val);
                }
            } else if (exist_B_1 && selectionRulesMomentum(state_row, state_col, 1) ) {
                real_t val = matrix_elements.getMagneticMomentum(state_row, state_col);
                if (std::abs(val) > tol) {
                    hamiltonian_magneticMomentum_p.addEntries(state_row.idx,state_col.idx,val);
                }
            } else if (exist_B_1 && selectionRulesMomentum(state_row, state_col, -1) ) {
                real_t val = matrix_elements.getMagneticMomentum(state_row, state_col);
                if (std::abs(val) > tol) {
                    hamiltonian_magneticMomentum_m.addEntries(state_row.idx,state_col.idx,val);
                }
            }

            if (diamagnetism && (exist_B_0 || exist_B_1) && selectionRulesMultipole(state_row, state_col, 0, 0)) {
                real_t val = matrix_elements.getDiamagnetism(state_row, state_col, 0);
                if (std::abs(val) > tol) {
                    hamiltonian_diamagnetism_00.addEntries(state_row.idx,state_col.idx,val);
                }
            } else if (diamagnetism && (exist_B_0 || exist_B_1) && selectionRulesMultipole(state_row, state_col, 2, 0)) {
                real_t val = matrix_elements.getDiamagnetism(state_row, state_col, 2);
                if (std::abs(val) > tol) {
                    hamiltonian_diamagnetism_20.addEntries(state_row.idx,state_col.idx,val);
                }
            } else if (diamagnetism && (exist_B_0 && exist_B_1) && selectionRulesMultipole(state_row, state_col, 2, 1)) {
                real_t val = matrix_elements.getDiamagnetism(state_row, state_col, 2);
                if (std::abs(val) > tol) {
                    hamiltonian_diamagnetism_2p.addEntries(state_row.idx,state_col.idx,val);
                }
            } else if (diamagnetism && (exist_B_0 && exist_B_1) && selectionRulesMultipole(state_row, state_col, 2, -1)) {
                real_t val = matrix_elements.getDiamagnetism(state_row, state_col, 2);
                if (std::abs(val) > tol) {
                    hamiltonian_diamagnetism_2m.addEntries(state_row.idx,state_col.idx,val);
                }
            } else if (diamagnetism && (exist_B_1) && selectionRulesMultipole(state_row, state_col, 2, 2)) {
                real_t val = matrix_elements.getDiamagnetism(state_row, state_col, 2);
                if (std::abs(val) > tol) {
                    hamiltonian_diamagnetism_2pp.addEntries(state_row.idx,state_col.idx,val);
                }
            } else if (diamagnetism && (exist_B_1) && selectionRulesMultipole(state_row, state_col, 2, -2)) {
                real_t val = matrix_elements.getDiamagnetism(state_row, state_col, 2);
                if (std::abs(val) > tol) {
                    hamiltonian_diamagnetism_2mm.addEntries(state_row.idx,state_col.idx,val);
                }
            }
        }
    }

    std::cout << "One-atom Hamiltonian, compress field Hamiltonian" << std::endl;

    hamiltonian_electricMomentum_0.compress(basis->dim(), basis->dim());
    hamiltonian_electricMomentum_p.compress(basis->dim(), basis->dim());
    hamiltonian_electricMomentum_m.compress(basis->dim(), basis->dim());

    hamiltonian_magneticMomentum_0.compress(basis->dim(), basis->dim());
    hamiltonian_magneticMomentum_p.compress(basis->dim(), basis->dim());
    hamiltonian_magneticMomentum_m.compress(basis->dim(), basis->dim());

    hamiltonian_diamagnetism_00.compress(basis->dim(), basis->dim());
    hamiltonian_diamagnetism_20.compress(basis->dim(), basis->dim());
    hamiltonian_diamagnetism_2p.compress(basis->dim(), basis->dim());
    hamiltonian_diamagnetism_2m.compress(basis->dim(), basis->dim());
    hamiltonian_diamagnetism_2pp.compress(basis->dim(), basis->dim());
    hamiltonian_diamagnetism_2mm.compress(basis->dim(), basis->dim());


    ////////////////////////////////////////////////////////
    ////// Prepare processing of Hamiltonians //////////////
    ////////////////////////////////////////////////////////

    // TODO Put the logic in its own class

    std::cout << "One-atom Hamiltonian, processe Hamiltonians" << std::endl;

    // === Open database ===
    boost::filesystem::path path_db;

    if (utils::is_complex<scalar_t>::value) {
        path_db = path_cache / "cache_matrix_complex.db";
    } else {
        path_db = path_cache / "cache_matrix_real.db";
    }
    sqlite::handle db(path_db.string());

    // === Initialize variables ===
    bool flag_perhapsmissingtable = true;

    matrix_path.resize(nSteps);
    matrix_diag.resize(nSteps); // TODO maybe remove
    params.resize(nSteps); // TODO maybe remove


    ////////////////////////////////////////////////////////
    ////// Loop through steps //////////////////////////////
    ////////////////////////////////////////////////////////

    publisher.send(">>TOT%7d", nSteps);

#pragma omp parallel for schedule(static, 1)

    // Loop through steps
    for (size_t step = 0; step < nSteps; ++step) {

        // === Get parameters for the current position inside the loop ===

        // Get fields
        real_t normalized_position = (nSteps > 1) ? step/(nSteps-1.) : 0;

        real_t Ex = min_E_x+normalized_position*(max_E_x-min_E_x);
        real_t Ey = min_E_y+normalized_position*(max_E_y-min_E_y);
        real_t Ez = min_E_z+normalized_position*(max_E_z-min_E_z);
        real_t Bx = min_B_x+normalized_position*(max_B_x-min_B_x);
        real_t By = min_B_y+normalized_position*(max_B_y-min_B_y);
        real_t Bz = min_B_z+normalized_position*(max_B_z-min_B_z);

        scalar_t E_0 = min_E_0+normalized_position*(max_E_0-min_E_0);
        scalar_t E_p = min_E_p+normalized_position*(max_E_p-min_E_p);
        scalar_t E_m = min_E_m+normalized_position*(max_E_m-min_E_m);
        scalar_t B_0 = min_B_0+normalized_position*(max_B_0-min_B_0);
        scalar_t B_p = min_B_p+normalized_position*(max_B_p-min_B_p);
        scalar_t B_m = min_B_m+normalized_position*(max_B_m-min_B_m);

        // Get configuration and save fields
        Configuration conf = basicconf;
        conf["Ex"] << Ex;
        conf["Ey"] << Ey;
        conf["Ez"] << Ez;
        conf["Bx"] << Bx;
        conf["By"] << By;
        conf["Bz"] << Bz;

        // === Create table if necessary ===
        std::stringstream query;
        std::string spacer = "";

        if (flag_perhapsmissingtable) {
            query << "CREATE TABLE IF NOT EXISTS cache_one (uuid text NOT NULL PRIMARY KEY, "
                     "created TIMESTAMP DEFAULT CURRENT_TIMESTAMP, "
                     "accessed TIMESTAMP DEFAULT CURRENT_TIMESTAMP";
            for (auto p: conf) {
                query << ", " << p.first << " text";
            }
            query << ", UNIQUE (";
            for (auto p: conf) {
                query << spacer << p.first;
                spacer = ", ";
            }
            query << "));";

            flag_perhapsmissingtable = false;
        }

        // === Get uuid as filename === // TODO put code in its own method
        std::string uuid = "";
        spacer = "";
        query << "SELECT uuid FROM cache_one WHERE ";
        for (auto p: conf) {
            query << spacer << p.first << "='" << p.second.str() << "'";
            spacer = " AND ";
        }
        query << ";";

#pragma omp critical(database)
        {
            sqlite::result result = db.query(query);
            if (result.size() == 1) {
                uuid = result.first();
            }
        }

        if (uuid != "") {
            query.str(std::string());
            query << "UPDATE cache_one SET accessed = CURRENT_TIMESTAMP WHERE uuid = '" << uuid << "';";
#pragma omp critical(database)
            db.exec(query.str()); // TODO check whether this slows down the program

        } else {
            boost::uuids::uuid u = generator();
            boost::algorithm::hex(u.begin(), u.end(), std::back_inserter(uuid));

            query.str(std::string());
            query << "INSERT INTO cache_one (uuid";
            for (auto p: conf) {
                query << ", " << p.first;
            }
            query << ") values ( '" << uuid << "'";
            for (auto p: conf) {
                query << ", " << "'" << p.second.str() << "'";
            }
            query << ");";
#pragma omp critical(database)
            db.exec(query);
        }

        // === Check existence of files === // TODO put code in its own method

        // Check whether .mat and .json file exists and compare settings in program with settings in .json file
        boost::filesystem::path path, path_mat, path_json;

        path = path_cache_mat / ("one_" + uuid);
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

        // calculate Hamiltonian if "is_existing" is false
        std::shared_ptr<Hamiltonianmatrix> mat;
        if (!is_existing || !totalmatrix.load(path_mat.string())) {

            // --- Build matrix ---
            totalmatrix = hamiltonian_energy
                    -hamiltonian_electricMomentum_0*E_0
                    +hamiltonian_electricMomentum_p*E_m
                    +hamiltonian_electricMomentum_m*E_p
                    +hamiltonian_magneticMomentum_0*B_0
                    -hamiltonian_magneticMomentum_p*B_m
                    -hamiltonian_magneticMomentum_m*B_p
                    +hamiltonian_diamagnetism_00*(B_0*B_0-2.*B_p*B_m)
                    -hamiltonian_diamagnetism_20*(B_0*B_0+B_p*B_m)
                    +std::sqrt(3)*hamiltonian_diamagnetism_2p*B_0*B_m
                    +std::sqrt(3)*hamiltonian_diamagnetism_2m*B_0*B_p
                    -std::sqrt(1.5)*hamiltonian_diamagnetism_2pp*B_m*B_m
                    -std::sqrt(1.5)*hamiltonian_diamagnetism_2mm*B_p*B_p;

            // Stdout: Hamiltonian assembled
#pragma omp critical(textoutput)
            {
                publisher.send(">>DIM%7d", totalmatrix.num_basisvectors());
                std::cout << "One-atom Hamiltonian, " <<  step+1 << ". Hamiltonian assembled" << std::endl;
            }

            // --- Diagonalize matrix and save diagonalized matrix ---
            totalmatrix.diagonalize();
            totalmatrix.save(path_mat.string());

            // Stdout: Hamiltonian diagonalized
#pragma omp critical(textoutput)
            {
                publisher.send(">>OUT%7d%7d%7d%7d %s", step+1, step, 1, 0, path.c_str());
                std::cout << "One-atom Hamiltonian, " <<  step+1 << ". Hamiltonian diagonalized" << std::endl;
            }
        } else {
            // Stdout: Hamiltonian loaded
#pragma omp critical(textoutput)
            {
                publisher.send(">>DIM%7d", totalmatrix.num_basisvectors());
                publisher.send(">>OUT%7d%7d%7d%7d %s", step+1, step, 1, 0, path.c_str());
                std::cout << "One-atom Hamiltonian, " <<  step+1 << ". Hamiltonian loaded" << std::endl;
            }
        }

        // === Store path to configuration and diagonalized matrix ===
        matrix_path[step] = path.string();
        matrix_diag[step] = std::make_shared<Hamiltonianmatrix>(totalmatrix); // TODO maybe remove
        params[step] = std::make_shared<Configuration>(conf); // TODO maybe remove
    }

    std::cout << "One-atom Hamiltonian, all Hamiltonians processed" << std::endl;

}
