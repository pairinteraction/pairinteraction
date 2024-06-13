#include "basis/BasisAtom.hpp"
#include "basis/BasisAtomCreator.hpp"
#include "database/Database.hpp"
#include "setup.hpp"
#include "system/SystemAtom.hpp"
#include "utils/args.hpp"

#include <filesystem>
#include <spdlog/spdlog.h>

int main(int argc, char **argv) {
    // Call the setup function to configure logging
    setup();

    // Create a database instance
    std::filesystem::path databasedir;
    bool download_missing = false;

    for (int i = 1; i < argc; ++i) {
        bool found = args::parse_download_missing(i, argc, argv, download_missing);
        if (!found) {
            found = args::parse_database(i, argc, argv, databasedir);
        }
    }

    Database database(download_missing, false, databasedir);

    // Create a basis
    auto basis = BasisAtomCreator<double>()
                     .set_species("Rb")
                     .restrict_quantum_number_n(60, 80)
                     .restrict_quantum_number_l(0, 10)
                     .create(database);

    SPDLOG_INFO("Number of basis states: {}", basis->get_number_of_states());

    // Create a system
    auto system = SystemAtom<double>(basis);
    system.set_electric_field({0, 0, 1});
    system.set_magnetic_field({0, 0, 1});
    system.enable_diamagnetism(true);

    return 0;
}
