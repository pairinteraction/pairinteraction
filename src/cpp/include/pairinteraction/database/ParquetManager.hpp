#pragma once

#include <filesystem>
#include <regex>
#include <string>
#include <unordered_map>
#include <vector>

namespace pairinteraction {

class GitHubDownloader;

class ParquetManager {
public:
    struct ParquetFileInfo {
        std::string path;
        int version = -1;
    };

    ParquetManager(std::filesystem::path directory, const GitHubDownloader &downloader,
                   std::vector<std::string> repo_paths);
    void scan_local();
    void scan_remote();
    std::string get_path(const std::string &table_name);
    std::string get_versions_info() const;

private:
    std::filesystem::path directory;
    const GitHubDownloader &downloader;
    std::vector<std::string> repo_paths;
    std::unordered_map<std::string, ParquetFileInfo> local_table_files;
    std::unordered_map<std::string, ParquetFileInfo> remote_table_files;
    std::regex file_regex = std::regex(R"(^(\w+)_v(\d+)\.parquet$)");
};

} // namespace pairinteraction
