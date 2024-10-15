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

    // Create a diagonalizer instance
    pairinteraction::DiagonalizerEigen<double> diagonalizer;

    // Create and diagonalize systems for two atoms
    auto basis1 = pairinteraction::BasisAtomCreator<double>()
                      .set_species("Rb")
                      .restrict_quantum_number_n(59, 61)
                      .restrict_quantum_number_l(0, 1)
                      .create(database);
    auto basis2 = pairinteraction::BasisAtomCreator<double>()
                      .set_species("Rb")
                      .restrict_quantum_number_n(59, 61)
                      .restrict_quantum_number_l(0, 1)
                      .create(database);

    pairinteraction::SystemAtom<double> system1(basis1);
    system1.set_electric_field({0, 0, 1.9446903811524456e-10});

    pairinteraction::SystemAtom<double> system2(basis2);
    system2.set_electric_field({0, 0, 1.9446903811524456e-10});

    pairinteraction::diagonalize<pairinteraction::SystemAtom<double>>({system1, system2},
                                                                      diagonalizer);

    // Create and diagonalize a combined systems around |60S; 60S>
    auto ket = pairinteraction::KetAtomCreator<double>()
                   .set_species("Rb")
                   .set_quantum_number_n(60)
                   .set_quantum_number_l(0)
                   .set_quantum_number_m(0.5)
                   .create(database);
    double min_energy = 2 * ket->get_energy() - 3 / 6579683.920501762;
    double max_energy = 2 * ket->get_energy() + 3 / 6579683.920501762;

    auto combined_basis = pairinteraction::BasisCombinedCreator<double>()
                              .add(system1)
                              .add(system2)
                              .restrict_energy(min_energy, max_energy)
                              .restrict_quantum_number_m(1, 1)
                              .create();

    std::vector<pairinteraction::SystemCombined<double>> combined_systems;
    combined_systems.reserve(5);
    for (int i = 1; i < 6; ++i) {
        pairinteraction::SystemCombined<double> system(combined_basis);
        system.set_interatomic_distance(i * 1e-6 / 5.29177210544e-11);
        combined_systems.push_back(std::move(system));
    }

    pairinteraction::diagonalize(combined_systems, diagonalizer);

    // Print information about the basis
    SPDLOG_INFO("Number of states in the basis of the first atom: {}",
                system1.get_basis()->get_number_of_states());
    SPDLOG_INFO("Number of states in the basis of the second atom: {}",
                system2.get_basis()->get_number_of_states());
    SPDLOG_INFO("Number of states in the combined basis: {}",
                combined_systems[0].get_basis()->get_number_of_states());

    std::stringstream ss;
    std::string separator = "";
    for (auto ket : *combined_basis) {
        ss << separator << *ket;
        separator = " | ";
    }
    SPDLOG_INFO("Kets in the combined basis: {}", ss.str());
}
