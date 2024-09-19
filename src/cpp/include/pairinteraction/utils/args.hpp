#pragma once

#include <filesystem>
#include <string>
#include <system_error>

namespace pairinteraction::args {
inline bool parse_download_missing(int &i, int argc, char **const argv, bool &download_missing) {
    if (i < argc && std::string(argv[i]) == "--download-missing") {
        download_missing = true;
        return true;
    }
    return false;
}

inline bool parse_database_dir(int &i, int argc, char **const argv,
                               std::filesystem::path &database_dir) {
    if (i < argc && std::string(argv[i]) == "--database-dir") {
        if (++i >= argc) {
            throw std::runtime_error(
                "Option --database-dir expects an arguments, but none is given.");
        }
        database_dir = argv[i];
        if (!std::filesystem::exists(database_dir)) {
            throw std::filesystem::filesystem_error(
                "Cannot access database", database_dir.string(),
                std::make_error_code(std::errc::no_such_file_or_directory));
        }
        database_dir = std::filesystem::canonical(database_dir);
        if (!std::filesystem::is_directory(database_dir)) {
            throw std::filesystem::filesystem_error(
                "Cannot access database", database_dir.string(),
                std::make_error_code(std::errc::not_a_directory));
        }
        return true;
    }
    return false;
}

inline bool parse_data_dir(int &i, int argc, char **const argv, std::filesystem::path &data_dir) {
    if (i < argc && std::string(argv[i]) == "--data-dir") {
        if (++i >= argc) {
            throw std::runtime_error("Option --data-dir expects an arguments, but none is given.");
        }
        data_dir = argv[i];
        if (!std::filesystem::exists(data_dir)) {
            throw std::filesystem::filesystem_error(
                "Cannot access data directory", data_dir.string(),
                std::make_error_code(std::errc::no_such_file_or_directory));
        }
        data_dir = std::filesystem::canonical(data_dir);
        if (!std::filesystem::is_directory(data_dir)) {
            throw std::filesystem::filesystem_error(
                "Cannot access data directory", data_dir.string(),
                std::make_error_code(std::errc::not_a_directory));
        }
        return true;
    }
    return false;
}

} // namespace pairinteraction::args
