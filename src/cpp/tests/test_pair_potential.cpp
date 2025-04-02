// SPDX-FileCopyrightText: 2024 Pairinteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "pairinteraction/pairinteraction.hpp"
#include "pairinteraction/utils/args.hpp"

#include <Eigen/Dense>
#include <filesystem>
#include <fstream>
#include <spdlog/spdlog.h>
#include <vector>

constexpr double HARTREE_IN_GHZ = 6579683.920501762;
constexpr double UM_IN_ATOMIC_UNITS = 1 / 5.29177210544e-5;

int main(int argc, char **argv) {
    // Call the setup function to configure logging
    pairinteraction::setup();

    // Create a database instance
    std::filesystem::path database_dir;
    std::filesystem::path data_dir;
    bool download_missing = false;

    for (int i = 1; i < argc; ++i) {
        bool found = pairinteraction::args::parse_download_missing(i, argc, argv, download_missing);
        if (!found) {
            found = pairinteraction::args::parse_database_dir(i, argc, argv, database_dir);
        }
        if (!found) {
            pairinteraction::args::parse_data_dir(i, argc, argv, data_dir);
        }
    }

    pairinteraction::Database database(download_missing, true, database_dir);

    // Create and diagonalize systems for two atoms
    auto basis = pairinteraction::BasisAtomCreator<double>()
                     .set_species("Rb")
                     .restrict_quantum_number_n(58, 62)
                     .restrict_quantum_number_l(0, 2)
                     .create(database);
    SPDLOG_INFO("Number of single-atom basis states: {}", basis->get_number_of_states());

    pairinteraction::SystemAtom<double> system(basis);

    // Create two-atom systems for different interatomic distances
    auto ket = pairinteraction::KetAtomCreator()
                   .set_species("Rb")
                   .set_quantum_number_n(60)
                   .set_quantum_number_l(0)
                   .set_quantum_number_m(0.5)
                   .create(database);
    double min_energy = 2 * ket->get_energy() - 3 / HARTREE_IN_GHZ;
    double max_energy = 2 * ket->get_energy() + 3 / HARTREE_IN_GHZ;

    auto basis_pair = pairinteraction::BasisPairCreator<double>()
                          .add(system)
                          .add(system)
                          .restrict_energy(min_energy, max_energy)
                          .restrict_quantum_number_m(1, 1)
                          .create();
    SPDLOG_INFO("Number of two-atom basis states: {}", basis_pair->get_number_of_states());

    std::vector<pairinteraction::SystemPair<double>> system_pairs;
    system_pairs.reserve(5);
    for (int i = 1; i < 6; ++i) {
        pairinteraction::SystemPair<double> system(basis_pair);
        system.set_distance_vector({0, 0, i * UM_IN_ATOMIC_UNITS});
        system_pairs.push_back(std::move(system));
    }

    // Diagonalize the systems in parallel
    pairinteraction::DiagonalizerEigen<double> diagonalizer(pairinteraction::FloatType::FLOAT32);
    pairinteraction::diagonalize(system_pairs, diagonalizer);

    // Sort by the eigenenergies
    for (auto &system : system_pairs) {
        auto sorter = system.get_sorter({pairinteraction::TransformationType::SORT_BY_ENERGY});
        system.transform(sorter);
    }

    // Extract results
    std::vector<std::string> kets;
    Eigen::MatrixX<double> eigenenergies(system_pairs.size(), basis_pair->get_number_of_states());
    Eigen::MatrixX<double> eigenstates(system_pairs.size(),
                                       basis_pair->get_number_of_states() *
                                           basis_pair->get_number_of_states());
    Eigen::MatrixX<double> overlaps(system_pairs.size(), basis_pair->get_number_of_states());

    kets.reserve(basis->get_number_of_states());
    for (const auto &ket : *basis_pair) {
        std::stringstream ss;
        ss << *ket;
        kets.push_back(ss.str());
    }

    for (Eigen::Index i = 0; i < static_cast<Eigen::Index>(system_pairs.size()); ++i) {
        eigenenergies.row(i) = system_pairs[i].get_eigenenergies() * HARTREE_IN_GHZ;

        Eigen::MatrixX<double> tmp =
            system_pairs[i].get_eigenbasis()->get_coefficients().toDense().transpose();
        eigenstates.row(i) = Eigen::Map<Eigen::VectorXd>(tmp.data(), tmp.size());

        overlaps.row(i) = system_pairs[i].get_eigenbasis()->get_overlaps(ket, ket);
    }

    // Compare with reference data
    bool success = true;

    // Check kets
    const std::filesystem::path reference_kets_file =
        data_dir / "reference_pair_potential/kets.txt";
    SPDLOG_INFO("Reference kets: {}", reference_kets_file.string());

    std::vector<std::string> reference_kets;
    std::string line;
    std::ifstream stream(reference_kets_file);
    if (!stream) {
        SPDLOG_ERROR("Could not open reference kets file");
        success = false;
    } else {
        reference_kets.reserve(basis_pair->get_number_of_states());
        while (std::getline(stream, line)) {
            reference_kets.push_back(line);
        }
        stream.close();
        if (kets.size() != reference_kets.size()) {
            SPDLOG_ERROR("Number of kets does not match reference data");
            success = false;
        } else if (kets != reference_kets) {
            for (size_t i = 0; i < kets.size(); ++i) {
                SPDLOG_DEBUG("Ket: {} vs {}, match: {}", kets[i], reference_kets[i],
                             kets[i] == reference_kets[i]);
            }
            SPDLOG_ERROR("Kets do not match reference data");
            success = false;
        }
    }

    // Check eigenenergies
    const std::filesystem::path reference_eigenenergies_file =
        data_dir / "reference_pair_potential/eigenenergies.txt";
    SPDLOG_INFO("Reference eigenenergies: {}", reference_eigenenergies_file.string());

    Eigen::MatrixXd reference_eigenenergies(system_pairs.size(),
                                            basis_pair->get_number_of_states());
    stream = std::ifstream(reference_eigenenergies_file);
    if (!stream) {
        SPDLOG_ERROR("Could not open reference eigenenergies file");
        success = false;
    } else {
        for (Eigen::Index i = 0; i < static_cast<Eigen::Index>(system_pairs.size()); ++i) {
            for (Eigen::Index j = 0;
                 j < static_cast<Eigen::Index>(basis_pair->get_number_of_states()); ++j) {
                stream >> reference_eigenenergies(i, j);
            }
        }
        stream.close();
        if ((eigenenergies - reference_eigenenergies).array().abs().maxCoeff() >
            100e-6) { // 100 kHz precision
            for (Eigen::Index i = 0; i < eigenenergies.size(); ++i) {
                SPDLOG_DEBUG("Eigenenergy: {} vs {}, delta: {}", eigenenergies(i),
                             reference_eigenenergies(i),
                             std::abs(eigenenergies(i) - reference_eigenenergies(i)));
            }
            SPDLOG_ERROR("Eigenenergies do not match reference data");
            success = false;
        }
    }

    // Check overlaps
    const std::filesystem::path reference_overlaps_file =
        data_dir / "reference_pair_potential/overlaps.txt";
    SPDLOG_INFO("Reference overlaps: {}", reference_overlaps_file.string());

    Eigen::MatrixXd reference_overlaps(system_pairs.size(), basis_pair->get_number_of_states());
    stream = std::ifstream(reference_overlaps_file);
    if (!stream) {
        SPDLOG_ERROR("Could not open reference overlap file");
        success = false;
    } else {
        for (Eigen::Index i = 0; i < static_cast<Eigen::Index>(system_pairs.size()); ++i) {
            for (Eigen::Index j = 0;
                 j < static_cast<Eigen::Index>(basis_pair->get_number_of_states()); ++j) {
                stream >> reference_overlaps(i, j);
            }
        }
        stream.close();
        if ((overlaps - reference_overlaps).array().abs().maxCoeff() > 1e-5) {
            for (Eigen::Index i = 0; i < overlaps.size(); ++i) {
                SPDLOG_DEBUG("Overlap: {} vs {}, delta: {}", overlaps(i), reference_overlaps(i),
                             std::abs(overlaps(i) - reference_overlaps(i)));
            }
            SPDLOG_ERROR("Overlaps do not match reference data");
            success = false;
        }
    }

    // Check eigenstates
    // Because of degeneracies, checking the eigenstates against reference data is complicated.
    // Thus, we only check their normalization and orthogonality.
    Eigen::VectorXd cumulative_norm =
        (eigenstates.adjoint().array() * eigenstates.transpose().array()).colwise().sum();
    if (!cumulative_norm.isApprox(Eigen::VectorXd::Constant(5, 19), 1e-5)) {
        SPDLOG_ERROR("Eigenvectors are not orthonormal.");
        success = false;
    }

    return success ? 0 : 1;
}
