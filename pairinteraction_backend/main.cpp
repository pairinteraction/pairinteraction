#include "setup.hpp"
#include "test.hpp"

#include <filesystem>
#include <string>
#include <system_error>
#include <vector>

int main(int argc, char **argv) {
    setup();

    std::filesystem::path databasedir = "";
    bool download_missing = false;
    std::vector<char *> args;
    args.reserve(argc);
    for (int i = 0; i < argc; ++i) {
        std::string key = argv[i];
        if (key == "--database") {
            if (++i >= argc) {
                throw std::runtime_error(
                    "Option --database expects an arguments, but none is given.");
            }
            databasedir = argv[i];
            if (!std::filesystem::exists(databasedir)) {
                throw std::filesystem::filesystem_error(
                    "Cannot access database", databasedir.string(),
                    std::make_error_code(std::errc::no_such_file_or_directory));
            }
            databasedir = std::filesystem::canonical(databasedir);
            if (!std::filesystem::is_directory(databasedir)) {
                throw std::filesystem::filesystem_error(
                    "Cannot access database", databasedir.string(),
                    std::make_error_code(std::errc::not_a_directory));
            }
        } else if (key == "--download-missing") {
            download_missing = true;
        } else {
            args.push_back(argv[i]);
        }
    }

    return test(args.size(), args.data(), download_missing, databasedir);
}
