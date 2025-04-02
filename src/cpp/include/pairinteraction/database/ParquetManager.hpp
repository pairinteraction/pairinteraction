// SPDX-FileCopyrightText: 2025 Pairinteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

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
    struct PathInfo {
        std::string path;
        bool cached = false;
    };

    struct LocalAssetInfo {
        int version_minor = -1;
        std::unordered_map<std::string, PathInfo> paths;
    };

    struct RemoteAssetInfo {
        int version_minor = -1;
        std::string endpoint;
    };

    ParquetManager(std::filesystem::path directory, const GitHubDownloader &downloader,
                   std::vector<std::string> repo_paths, duckdb::Connection &con, bool use_cache);
    void scan_local();
    void scan_remote();
    std::string get_path(const std::string &key, const std::string &table);
    std::string get_versions_info() const;

private:
    void react_on_rate_limit_reached(std::time_t reset_time);
    void update_local_asset(const std::string &key);
    void cache_table(std::unordered_map<std::string, PathInfo>::iterator table_it);

    std::filesystem::path directory_;
    const GitHubDownloader &downloader;
    std::vector<std::string> repo_paths_;
    duckdb::Connection &con;
    bool use_cache_;
    std::unordered_map<std::string, LocalAssetInfo> local_asset_info;
    std::unordered_map<std::string, RemoteAssetInfo> remote_asset_info;
    std::regex local_regex{R"(^(\w+)_v(\d+)\.(\d+)$)"};
    std::regex remote_regex{R"(^(\w+)_v(\d+)\.(\d+)\.zip$)"};
    std::shared_mutex mtx_local;
};

} // namespace pairinteraction
