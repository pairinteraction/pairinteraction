#include "pintr/pintr.hpp"
#include "pintr/utils/args.hpp"

#include <filesystem>
#include <spdlog/spdlog.h>
#include <vector>

int main(int argc, char **argv) {
    // Call the setup function to configure logging
    pintr::setup();

    // Create a database instance
    std::filesystem::path databasedir;
    bool download_missing = false;

    for (int i = 1; i < argc; ++i) {
        bool found = pintr::args::parse_download_missing(i, argc, argv, download_missing);
        if (!found) {
            pintr::args::parse_database(i, argc, argv, databasedir);
        }
    }

    thread_local pintr::Database database(download_missing, true, databasedir);

    // Create a basis
    auto basis = pintr::BasisAtomCreator<double>()
                     .set_species("Rb")
                     .restrict_quantum_number_n(58, 62)
                     .restrict_quantum_number_l(0, 2)
                     .restrict_quantum_number_m(0.5, 0.5)
                     .create(database);

    SPDLOG_INFO("Number of basis states: {}", basis->get_number_of_states());

    // Create systems for different values of the electric field
    std::vector<pintr::SystemAtom<double>> systems;
    systems.reserve(10);
    for (int i = 0; i < 10; ++i) {
        auto system = pintr::SystemAtom<double>(basis);
        system.set_electric_field({0, 0, i * 1e-9});
        systems.push_back(std::move(system));
    }

    // Diagonalize the systems in parallel
    pintr::DiagonalizerEigen<double> diagonalizer;
    pintr::diagonalize(systems, diagonalizer);

    return 0;
}
