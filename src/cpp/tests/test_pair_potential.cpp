#include "pairinteraction/pairinteraction.hpp"
#include "pairinteraction/utils/args.hpp"

#include <Eigen/Dense>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <spdlog/spdlog.h>
#include <vector>

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

    thread_local pairinteraction::Database database(download_missing, true, database_dir);

    // Create and diagonalize systems for two atoms
    auto basis = pairinteraction::BasisAtomCreator<double>()
                     .set_species("Rb")
                     .restrict_quantum_number_n(58, 62)
                     .restrict_quantum_number_l(0, 2)
                     .create(database);
    SPDLOG_INFO("Number of single-atom basis states: {}", basis->get_number_of_states());

    pairinteraction::SystemAtom<double> system(basis);

    // Create two-atom systems for different interatomic distances
    auto ket = pairinteraction::KetAtomCreator<double>()
                   .set_species("Rb")
                   .set_quantum_number_n(60)
                   .set_quantum_number_l(0)
                   .set_quantum_number_m(0.5)
                   .create(database);
    double min_energy = 2 * ket->get_energy() - 3 / 6579683.920501762;
    double max_energy = 2 * ket->get_energy() + 3 / 6579683.920501762;

    auto combined_basis = pairinteraction::BasisCombinedCreator<double>()
                              .add(system)
                              .add(system)
                              .restrict_energy(min_energy, max_energy)
                              .restrict_quantum_number_m(1, 1)
                              .create();
    SPDLOG_INFO("Number of two-atom basis states: {}", combined_basis->get_number_of_states());

    std::vector<pairinteraction::SystemCombined<double>> combined_systems;
    combined_systems.reserve(5);
    for (int i = 1; i < 6; ++i) {
        pairinteraction::SystemCombined<double> system(combined_basis);
        system.set_distance(i * 1e-6 / 5.29177210544e-11);
        combined_systems.push_back(std::move(system));
    }

    // Diagonalize the systems in parallel
    pairinteraction::DiagonalizerEigen<double> diagonalizer;
    pairinteraction::diagonalize(combined_systems, diagonalizer);

    // Sort by the eigenvalues
    for (auto &system : combined_systems) {
        auto sorter = system.get_sorter({pairinteraction::TransformationType::SORT_BY_ENERGY});
        system.transform(sorter);
    }

    // Extract results
    std::vector<std::string> kets;
    Eigen::MatrixX<double> eigenvalues(combined_systems.size(),
                                       combined_basis->get_number_of_states());
    Eigen::MatrixX<double> eigenstates(combined_systems.size(),
                                       combined_basis->get_number_of_states() *
                                           combined_basis->get_number_of_states());
    Eigen::MatrixX<double> overlaps(combined_systems.size(),
                                    combined_basis->get_number_of_states());

    kets.reserve(basis->get_number_of_states());
    for (const auto &ket : *combined_basis) {
        std::stringstream ss;
        ss << *ket;
        kets.push_back(ss.str());
    }

    for (Eigen::Index i = 0; i < static_cast<Eigen::Index>(combined_systems.size()); ++i) {
        eigenvalues.row(i) = combined_systems[i].get_eigenvalues() * 6579683.920501762;

        Eigen::MatrixX<double> tmp =
            combined_systems[i].get_eigenbasis()->get_coefficients().toDense().transpose();
        eigenstates.row(i) = Eigen::Map<Eigen::VectorXd>(tmp.data(), tmp.size());

        overlaps.row(i) = combined_systems[i].get_eigenbasis()->get_overlaps(ket, ket);
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
        reference_kets.reserve(combined_basis->get_number_of_states());
        while (std::getline(stream, line)) {
            reference_kets.push_back(line);
        }
        stream.close();
        if (kets != reference_kets) {
            SPDLOG_ERROR("Kets do not match reference data");
            success = false;
        }
    }

    // Check eigenvalues
    const std::filesystem::path reference_eigenvalues_file =
        data_dir / "reference_pair_potential/eigenvalues.txt";
    SPDLOG_INFO("Reference eigenvalues: {}", reference_eigenvalues_file.string());

    Eigen::MatrixXd reference_eigenvalues(combined_systems.size(),
                                          combined_basis->get_number_of_states());
    stream = std::ifstream(reference_eigenvalues_file);
    if (!stream) {
        SPDLOG_ERROR("Could not open reference eigenvalues file");
        success = false;
    } else {
        for (Eigen::Index i = 0; i < static_cast<Eigen::Index>(combined_systems.size()); ++i) {
            for (Eigen::Index j = 0;
                 j < static_cast<Eigen::Index>(combined_basis->get_number_of_states()); ++j) {
                stream >> reference_eigenvalues(i, j);
            }
        }
        stream.close();
        if (!eigenvalues.isApprox(reference_eigenvalues, 1e-9)) {
            for (Eigen::Index i = 0; i < eigenvalues.size(); ++i) {
                SPDLOG_DEBUG("Eigenvalue: {} vs {}, delta: {}", eigenvalues(i),
                             reference_eigenvalues(i),
                             std::abs(eigenvalues(i) - reference_eigenvalues(i)));
            }
            SPDLOG_ERROR("Eigenvalues do not match reference data");
            success = false;
        }
    }

    // Check overlaps
    const std::filesystem::path reference_overlaps_file =
        data_dir / "reference_pair_potential/overlaps.txt";
    SPDLOG_INFO("Reference overlaps: {}", reference_overlaps_file.string());

    Eigen::MatrixXd reference_overlaps(combined_systems.size(),
                                       combined_basis->get_number_of_states());
    stream = std::ifstream(reference_overlaps_file);
    if (!stream) {
        SPDLOG_ERROR("Could not open reference overlap file");
        success = false;
    } else {
        for (Eigen::Index i = 0; i < static_cast<Eigen::Index>(combined_systems.size()); ++i) {
            for (Eigen::Index j = 0;
                 j < static_cast<Eigen::Index>(combined_basis->get_number_of_states()); ++j) {
                stream >> reference_overlaps(i, j);
            }
        }
        stream.close();
        if (!overlaps.isApprox(reference_overlaps, 1e-9)) {
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
    if (!cumulative_norm.isApprox(Eigen::VectorXd::Constant(5, 19))) {
        SPDLOG_ERROR("Eigenvectors are not orthonormal.");
        success = false;
    }

    return success ? 0 : 1;
}
