#pragma once

#include <filesystem>
#include <regex>
#include <shared_mutex>
#include <string>
#include <unordered_map>
#include <vector>

namespace duckdb {
class Connection;
} // namespace duckdb

namespace pairinteraction {

class GitHubDownloader;

class ParquetManager {
public:
    struct LocalPathInfo {
        std::string path;
        bool cached = false;
    };

    ParquetManager(std::filesystem::path directory, const GitHubDownloader &downloader,
                   std::vector<std::string> repo_paths, duckdb::Connection &con, bool use_cache);
    void scan_local();
    void scan_remote();
    std::string get_path(const std::string &species, const std::string &table);
    std::string get_versions_info() const;

private:
    void react_on_rate_limit_reached(std::time_t reset_time);
    std::unordered_map<std::string, LocalPathInfo>::iterator
    update_table(const std::string &species, const std::string &table);
    void cache_table(std::unordered_map<std::string, LocalPathInfo>::iterator local_it);

    std::filesystem::path directory_;
    const GitHubDownloader &downloader;
    std::vector<std::string> repo_paths_;
    duckdb::Connection &con;
    bool use_cache_;
    std::unordered_map<std::string, int> table_local_versions;
    std::unordered_map<std::string, int> table_remote_versions;
    std::unordered_map<std::string, LocalPathInfo> table_local_paths;
    std::unordered_map<std::string, std::string> table_remote_paths; // TODO combine
    std::regex web_regex = std::regex(R"(^(\w+)_v(\d+)\.(\d+)\.(?:parquet|zip)$)");
    std::regex file_regex = std::regex(R"(^(\w+)_v(\d+)\.(\d+)\.parquet$)");
    std::regex dir_regex = std::regex(R"(^(\w+)_v(\d+)\.(\d+)$)");
    std::shared_mutex mtx_local_table_files;
};

} // namespace pairinteraction
