#include "pairinteraction/tools/setup.hpp"
#include "pairinteraction/tools/test.hpp"
#include "pairinteraction/utils/args.hpp"

#include <filesystem>
#include <vector>

int main(int argc, char **argv) {
    pairinteraction::setup();

    std::filesystem::path database_dir;
    bool download_missing = false;
    std::vector<char *> args;
    args.reserve(argc);

    for (int i = 0; i < argc; ++i) {
        bool found = pairinteraction::args::parse_download_missing(i, argc, argv, download_missing);
        if (!found) {
            found = pairinteraction::args::parse_database_dir(i, argc, argv, database_dir);
        }
        if (!found) {
            args.push_back(argv[i]);
        }
    }

    return pairinteraction::test(static_cast<int>(args.size()), args.data(), download_missing,
                                 database_dir);
}
