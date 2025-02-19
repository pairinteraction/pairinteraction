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
        int version;
    };

    ParquetManager(const std::filesystem::path &directory, GitHubDownloader &downloader,
                   const std::vector<std::string> &repo_paths);
    void scan_local();
    void scan_remote();
    std::string get_path(const std::string &table_name);

private:
    std::filesystem::path directory_;
    GitHubDownloader &downloader_;
    std::vector<std::string> repo_paths_;
    std::unordered_map<std::string, ParquetFileInfo> local_table_files_;
    std::unordered_map<std::string, ParquetFileInfo> remote_table_files_;
    std::regex file_regex_ = std::regex(R"(^(\w+)_v(\d+)\.parquet$)");
};

} // namespace pairinteraction
