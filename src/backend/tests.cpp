#include "pintr/tools/setup.hpp"
#include "pintr/tools/test.hpp"
#include "pintr/utils/args.hpp"

#include <filesystem>
#include <vector>

int main(int argc, char **argv) {
    setup();

    std::filesystem::path databasedir;
    bool download_missing = false;
    std::vector<char *> args;
    args.reserve(argc);

    for (int i = 0; i < argc; ++i) {
        bool found = args::parse_download_missing(i, argc, argv, download_missing);
        if (!found) {
            found = args::parse_database(i, argc, argv, databasedir);
        }
        if (!found) {
            args.push_back(argv[i]);
        }
    }

    return test(static_cast<int>(args.size()), args.data(), download_missing, databasedir);
}
