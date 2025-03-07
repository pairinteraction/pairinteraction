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
    struct AssetInfo {
        std::string path;
        bool is_dir = false;
        int version_minor = -1;
    };

    struct PathInfo {
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
    std::unordered_map<std::string, PathInfo>::iterator update_table(const std::string &species,
                                                                     const std::string &table);
    void cache_table(std::unordered_map<std::string, PathInfo>::iterator path_it);

    std::filesystem::path directory_;
    const GitHubDownloader &downloader;
    std::vector<std::string> repo_paths_;
    duckdb::Connection &con;
    bool use_cache_;
    std::unordered_map<std::string, AssetInfo> remote_asset_info;
    std::unordered_map<std::string, AssetInfo> local_asset_info;
    std::unordered_map<std::string, PathInfo> local_table_info;
    std::regex parquet_regex{R"(^(\w+)_v(\d+)\.(\d+)\.parquet$)"};
    std::regex zip_regex{R"(^(\w+)_v(\d+)\.(\d+)\.zip$)"};
    std::regex dir_regex{R"(^(\w+)_v(\d+)\.(\d+)$)"};
    std::shared_mutex mtx_local;
};

} // namespace pairinteraction
