#pragma once

#include <filesystem>
#include <string>
#include <system_error>

namespace args {
inline bool parse_download_missing(int &i, int argc, char **const argv, bool &download_missing) {
    if (i < argc && std::string(argv[i]) == "--download-missing") {
        download_missing = true;
        return true;
    }
    return false;
}

inline bool parse_database(int &i, int argc, char **const argv,
                           std::filesystem::path &databasedir) {
    if (i < argc && std::string(argv[i]) == "--database") {
        if (++i >= argc) {
            throw std::runtime_error("Option --database expects an arguments, but none is given.");
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
        return true;
    }
    return false;
}

} // namespace args
