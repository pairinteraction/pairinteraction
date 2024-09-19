#include "pairinteraction/pairinteraction.hpp"
#include "pairinteraction/utils/args.hpp"

#include <filesystem>
#include <spdlog/spdlog.h>

int main(int argc, char **argv) {
    // Call the setup function to configure logging
    pairinteraction::setup();

    // Create a database instance
    std::filesystem::path database_dir;
    bool download_missing = false;

    for (int i = 1; i < argc; ++i) {
        bool found = pairinteraction::args::parse_download_missing(i, argc, argv, download_missing);
        if (!found) {
            pairinteraction::args::parse_database(i, argc, argv, database_dir);
        }
    }

    pairinteraction::Database database(download_missing, true, database_dir);

    // Create a basis
    auto basis = pairinteraction::BasisAtomCreator<double>()
                     .set_species("Rb")
                     .restrict_quantum_number_n(58, 62)
                     .restrict_quantum_number_l(0, 3)
                     .create(database);

    SPDLOG_INFO("Number of basis states: {}", basis->get_number_of_states());

    // Create a system
    auto system = pairinteraction::SystemAtom<double>(basis);
    system.set_electric_field({0, 0, 1});
    system.set_magnetic_field({0, 0, 1});
    system.enable_diamagnetism(true);

    return 0;
}
