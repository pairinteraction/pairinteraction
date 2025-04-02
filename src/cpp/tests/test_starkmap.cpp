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
constexpr double VOLT_PER_CM_IN_ATOMIC_UNITS = 1 / 5.14220675112e9;

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

    // Create a basis
    auto ket = pairinteraction::KetAtomCreator()
                   .set_species("Rb")
                   .set_quantum_number_n(60)
                   .set_quantum_number_l(0)
                   .set_quantum_number_m(0.5)
                   .create(database);

    auto basis = pairinteraction::BasisAtomCreator<double>()
                     .set_species("Rb")
                     .restrict_quantum_number_n(58, 62)
                     .restrict_quantum_number_l(0, 2)
                     .create(database);

    SPDLOG_INFO("Number of basis states: {}", basis->get_number_of_states());

    // Create systems for different values of the electric field
    std::vector<pairinteraction::SystemAtom<double>> systems;
    systems.reserve(11);
    for (int i = 0; i < 11; ++i) {
        pairinteraction::SystemAtom<double> system(basis);
        system.set_electric_field({0, 0, i * VOLT_PER_CM_IN_ATOMIC_UNITS});
        systems.push_back(std::move(system));
    }

    // Diagonalize the systems in parallel
    pairinteraction::DiagonalizerEigen<double> diagonalizer(pairinteraction::FloatType::FLOAT32);
    pairinteraction::diagonalize(systems, diagonalizer);

    // Sort by the eigenenergies
    for (auto &system : systems) {
        auto sorter = system.get_sorter({pairinteraction::TransformationType::SORT_BY_ENERGY});
        system.transform(sorter);
    }

    // Extract results
    std::vector<std::string> kets;
    Eigen::MatrixX<double> eigenenergies(systems.size(), basis->get_number_of_states());
    Eigen::MatrixX<double> eigenstates(
        systems.size(), basis->get_number_of_states() * basis->get_number_of_states());
    Eigen::MatrixX<double> overlaps(systems.size(), basis->get_number_of_states());

    kets.reserve(basis->get_number_of_states());
    for (const auto &ket : *basis) {
        std::stringstream ss;
        ss << *ket;
        kets.push_back(ss.str());
    }

    for (Eigen::Index i = 0; i < static_cast<Eigen::Index>(systems.size()); ++i) {
        eigenenergies.row(i) = systems[i].get_eigenenergies() * HARTREE_IN_GHZ;

        Eigen::MatrixX<double> tmp =
            systems[i].get_eigenbasis()->get_coefficients().toDense().transpose();
        eigenstates.row(i) = Eigen::Map<Eigen::VectorXd>(tmp.data(), tmp.size());

        overlaps.row(i) = systems[i].get_eigenbasis()->get_overlaps(ket);
    }

    // Compare with reference data
    bool success = true;

    // Check kets
    const std::filesystem::path reference_kets_file = data_dir / "reference_stark_map/kets.txt";
    SPDLOG_INFO("Reference kets: {}", reference_kets_file.string());

    std::vector<std::string> reference_kets;
    std::string line;
    std::ifstream stream(reference_kets_file);
    if (!stream) {
        SPDLOG_ERROR("Could not open reference kets file");
        success = false;
    } else {
        reference_kets.reserve(basis->get_number_of_states());
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
        data_dir / "reference_stark_map/eigenenergies.txt";
    SPDLOG_INFO("Reference eigenenergies: {}", reference_eigenenergies_file.string());

    Eigen::MatrixXd reference_eigenenergies(systems.size(), basis->get_number_of_states());
    stream = std::ifstream(reference_eigenenergies_file);
    if (!stream) {
        SPDLOG_ERROR("Could not open reference eigenenergies file");
        success = false;
    } else {
        for (Eigen::Index i = 0; i < static_cast<Eigen::Index>(systems.size()); ++i) {
            for (Eigen::Index j = 0; j < static_cast<Eigen::Index>(basis->get_number_of_states());
                 ++j) {
                stream >> reference_eigenenergies(i, j);
            }
        }
        stream.close();
        if ((eigenenergies - reference_eigenenergies).array().abs().maxCoeff() >
            500e-6) { // 500 kHz precision
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
        data_dir / "reference_stark_map/overlaps.txt";
    SPDLOG_INFO("Reference overlaps: {}", reference_overlaps_file.string());

    Eigen::MatrixXd reference_overlaps(systems.size(), basis->get_number_of_states());
    stream = std::ifstream(reference_overlaps_file);
    if (!stream) {
        SPDLOG_ERROR("Could not open reference overlap file");
        success = false;
    } else {
        for (Eigen::Index i = 0; i < static_cast<Eigen::Index>(systems.size()); ++i) {
            for (Eigen::Index j = 0; j < static_cast<Eigen::Index>(basis->get_number_of_states());
                 ++j) {
                stream >> reference_overlaps(i, j);
            }
        }
        stream.close();
        if ((overlaps - reference_overlaps).array().abs().maxCoeff() > 5e-5) {
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
    if (!cumulative_norm.isApprox(Eigen::VectorXd::Constant(11, 90), 5e-5)) {
        SPDLOG_ERROR("Eigenvectors are not orthonormal.");
        success = false;
    }

    return success ? 0 : 1;
}
