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
    struct ParquetFileInfo {
        std::string path;
        int version = -1;
        bool cached = false;
    };

    ParquetManager(std::filesystem::path directory, const GitHubDownloader &downloader,
                   std::vector<std::string> repo_paths, duckdb::Connection &con, bool use_cache);
    void scan_local();
    void scan_remote();
    std::string get_path(const std::string &table_name);
    std::string get_versions_info() const;

private:
    void react_on_rate_limit_reached(std::time_t reset_time);
    std::unordered_map<std::string, ParquetFileInfo>::iterator
    update_table(const std::string &table_name);
    void cache_table(std::unordered_map<std::string, ParquetFileInfo>::iterator local_it);

    std::filesystem::path directory_;
    const GitHubDownloader &downloader;
    std::vector<std::string> repo_paths_;
    duckdb::Connection &con;
    bool use_cache_;
    std::unordered_map<std::string, ParquetFileInfo> local_table_files;
    std::unordered_map<std::string, ParquetFileInfo> remote_table_files;
    std::regex file_regex = std::regex(R"(^(\w+)_v(\d+)\.parquet$)");
    std::shared_mutex mtx_local_table_files;
};

} // namespace pairinteraction
