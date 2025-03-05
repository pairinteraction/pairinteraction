#include "pairinteraction/database/ParquetManager.hpp"

#include "pairinteraction/database/GitHubDownloader.hpp"

#include <cpptrace/cpptrace.hpp>
#include <ctime>
#include <duckdb.hpp>
#include <filesystem>
#include <fmt/core.h>
#include <fstream>
#include <future>
#include <iomanip>
#include <nlohmann/json.hpp>
#include <regex>
#include <set>
#include <spdlog/spdlog.h>
#include <sstream>
#include <stdexcept>
#include <string>

std::string format_time(std::time_t time_val) {
    std::tm *ptm = std::localtime(&time_val);
    std::ostringstream oss;
    oss << std::put_time(ptm, "%Y-%m-%d %H:%M:%S");
    return oss.str();
}

namespace pairinteraction {
ParquetManager::ParquetManager(std::filesystem::path directory, const GitHubDownloader &downloader,
                               std::vector<std::string> repo_paths, duckdb::Connection &con,
                               bool use_cache)
    : directory_(std::move(directory)), downloader(downloader), repo_paths_(std::move(repo_paths)),
      con(con), use_cache_(use_cache) {
    // Ensure the local directory exists
    if (!std::filesystem::exists(directory_)) {
        std::filesystem::create_directories(directory_);
    }

    // If repo paths are provided, check the GitHub rate limit
    if (!repo_paths_.empty()) {
        auto rate_limit = downloader.get_rate_limit();
        if (rate_limit.remaining <= 0) {
            react_on_rate_limit_reached(rate_limit.reset_time);
        } else {
            SPDLOG_INFO("Remaining GitHub API requests: {}. Rate limit resets at {}.",
                        rate_limit.remaining, format_time(rate_limit.reset_time));
        }
    }
}

void ParquetManager::scan_remote() {
    remote_table_files.clear();

    std::vector<std::future<GitHubDownloader::Result>> futures;
    std::vector<std::filesystem::path> cache_files;
    std::vector<nlohmann::json> cached_docs;

    futures.reserve(repo_paths_.size());
    cache_files.reserve(repo_paths_.size());
    cached_docs.reserve(repo_paths_.size());

    // For each repo path, load its cached JSON (or an empty JSON) and issue the download
    for (const auto &remote_endpoint : repo_paths_) {
        // Generate a unique cache filename per endpoint
        std::string cache_filename =
            "homepage_cache_" + std::to_string(std::hash<std::string>{}(remote_endpoint)) + ".json";
        std::filesystem::path cache_file = directory_ / cache_filename;
        cache_files.push_back(cache_file);

        // Load cached JSON from file if it exists
        nlohmann::json cached_doc;
        if (std::filesystem::exists(cache_file)) {
            try {
                std::ifstream in(cache_file);
                in.exceptions(std::ifstream::failbit | std::ifstream::badbit);
                in >> cached_doc;
            } catch (const std::exception &e) {
                SPDLOG_WARN("Error reading or parsing {}: {}; using empty JSON.",
                            cache_file.string(), e.what());
                cached_doc = nlohmann::json{};
            }
        }
        cached_docs.push_back(cached_doc);

        // Extract the last-modified header from the cached JSON if available
        std::string last_modified;
        if (!cached_doc.is_null() && cached_doc.contains("last-modified")) {
            last_modified = cached_doc["last-modified"].get<std::string>();
        }

        // Issue the asynchronous download using the cached last-modified value.
        futures.push_back(downloader.download(remote_endpoint, last_modified));
    }

    // Process downloads for each repo path
    for (size_t i = 0; i < futures.size(); ++i) {
        const auto &remote_endpoint = repo_paths_[i];
        auto result = futures[i].get();
        const auto &cache_file = cache_files[i];
        const auto &cached_doc = cached_docs[i];

        nlohmann::json doc;
        if (result.status_code == 200) {
            doc = nlohmann::json::parse(result.body, nullptr, /*allow_exceptions=*/false);
            doc["last-modified"] = result.last_modified;
            std::ofstream out(cache_file);
            if (!out) {
                throw std::runtime_error(
                    fmt::format("Failed to open {} for writing", cache_file.string()));
            }
            out << doc;
            out.close();
            SPDLOG_INFO("Using downloaded overview of available tables from {}.", remote_endpoint);
        } else if (result.status_code == 304) {
            if (cached_doc.is_null() || cached_doc.empty()) {
                throw std::runtime_error(
                    fmt::format("Received 304 Not Modified but cached response {} does not exist.",
                                cache_file.string()));
            }
            doc = cached_doc;
            SPDLOG_INFO("Using cached overview of available tables from {}.", remote_endpoint);
        } else if (result.status_code == 403 || result.status_code == 429) {
            react_on_rate_limit_reached(result.rate_limit.reset_time);
            return;
        } else {
            throw std::runtime_error(fmt::format(
                "Failed to download overview of available tables from {}: status code {}.",
                remote_endpoint, result.status_code));
        }

        // Validate the JSON response
        if (doc.is_discarded() || !doc.contains("assets")) {
            throw std::runtime_error(fmt::format(
                "Failed to parse remote JSON or missing 'assets' key from {}.", remote_endpoint));
        }

        // Update remote_table_files based on the asset entries
        for (auto &asset : doc["assets"]) {
            std::string filename = asset["name"].get<std::string>();
            std::smatch match;
            if (std::regex_match(filename, match, file_regex) && match.size() == 3) {
                std::string table_name = match[1].str();
                int remote_version = std::stoi(match[2].str());
                auto it = remote_table_files.find(table_name);
                if (it == remote_table_files.end() || remote_version > it->second.version) {
                    std::string remote_url = asset["url"].get<std::string>();
                    const std::string host = downloader.get_host();
                    remote_table_files[table_name] = {remote_url.erase(0, host.size()),
                                                      remote_version, false};
                }
            }
        }
    }
}

void ParquetManager::scan_local() {
    local_table_files.clear();

    // Iterate over files in the directory to update local_table_files
    for (const auto &entry : std::filesystem::directory_iterator(directory_)) {
        if (entry.is_regular_file()) {
            std::string filename = entry.path().filename().string();
            std::smatch match;
            if (std::regex_match(filename, match, file_regex) && match.size() == 3) {
                std::string table_name = match[1].str();
                int version = std::stoi(match[2].str());
                auto it = local_table_files.find(table_name);
                if (it == local_table_files.end() || version > it->second.version) {
                    local_table_files[table_name] = {entry.path().string(), version, false};
                }
            }
        }
    }
}

void ParquetManager::react_on_rate_limit_reached(std::time_t reset_time) {
    repo_paths_.clear();
    remote_table_files.clear();
    SPDLOG_WARN("Rate limit reached, resets at {}. The download of database tables is disabled.",
                format_time(reset_time));
}

std::unordered_map<std::string, ParquetManager::ParquetFileInfo>::iterator
ParquetManager::update_table(const std::string &table_name) {
    // Get remote version if available
    int remote_version = -1;
    auto remote_it = remote_table_files.find(table_name);
    if (remote_it != remote_table_files.end()) {
        remote_version = remote_it->second.version;
    }

    // Get local version if available and check if it is up-to-date
    {
        int local_version = -1;
        std::shared_lock<std::shared_mutex> lock(mtx_local_table_files);
        auto local_it = local_table_files.find(table_name);
        if (local_it != local_table_files.end()) {
            local_version = local_it->second.version;
        }
        if (local_version >= remote_version) {
            return local_it;
        }
    }

    // If it is not up-to-date, acquire a unique lock for updating the table
    std::unique_lock<std::shared_mutex> lock(mtx_local_table_files);

    // Re-check if the table is up to date because another thread might have updated it
    int local_version = -1;
    auto local_it = local_table_files.find(table_name);
    if (local_it != local_table_files.end()) {
        local_version = local_it->second.version;
    }
    if (local_version >= remote_version) {
        return local_it;
    }

    SPDLOG_INFO("Updating table {} from version {} to version {}.", table_name, local_version,
                remote_version);

    // Download the remote file
    auto result = downloader.download(remote_it->second.path, "", true).get();
    if (result.status_code == 403 || result.status_code == 429) {
        react_on_rate_limit_reached(result.rate_limit.reset_time);
        return local_it;
    }
    if (result.status_code != 200) {
        throw std::runtime_error(fmt::format("Failed to download table {}: status code {}.",
                                             table_name, result.status_code));
    }

    // Save the downloaded file locally
    std::string local_file_path =
        (directory_ / (table_name + "_v" + std::to_string(remote_version) + ".parquet")).string();
    std::ofstream out(local_file_path, std::ios::binary);
    if (!out) {
        throw std::runtime_error(fmt::format("Failed to open {} for writing", local_file_path));
    }
    out << result.body;
    out.close();

    // Update the internal mapping
    return local_table_files
        .insert_or_assign(table_name, ParquetFileInfo{local_file_path, remote_version, false})
        .first;
}

void ParquetManager::cache_table(
    std::unordered_map<std::string, ParquetFileInfo>::iterator local_it) {
    // Check if the table is already cached
    {
        std::shared_lock<std::shared_mutex> lock(mtx_local_table_files);
        if (local_it->second.cached) {
            return;
        }
    }

    // Acquire a unique lock for caching the table
    std::unique_lock<std::shared_mutex> lock(mtx_local_table_files);

    // Re-check if the table is already cached because another thread might have cached it
    if (local_it->second.cached) {
        return;
    }

    // Cache the table in memory
    auto result = con.Query(fmt::format(R"(CREATE TEMP TABLE '{}' AS SELECT * FROM '{}')",
                                        local_it->first, local_it->second.path));
    if (result->HasError()) {
        throw cpptrace::runtime_error("Error creating table: " + result->GetError());
    }
    local_it->second.path = local_it->first;
    local_it->second.cached = true;
}

std::string ParquetManager::get_path(const std::string &table_name) {
    // Update the local table if a newer version is available remotely
    auto local_it = this->update_table(table_name);

    // Cache the local table in memory if requested
    if (use_cache_) {
        this->cache_table(local_it);
    }

    // Return the path to the local table file
    if (local_it != local_table_files.end()) {
        return local_it->second.path;
    }
    throw std::runtime_error("Table " + table_name + " not found.");
}

std::string ParquetManager::get_versions_info() const {
    // Helper lambda returns the version string if available
    auto get_version = [](const auto &map, const std::string &table) -> int {
        if (auto it = map.find(table); it != map.end()) {
            return it->second.version;
        }
        return -1;
    };

    // Gather all unique table names
    std::set<std::string> tables;
    for (const auto &entry : local_table_files) {
        tables.insert(entry.first);
    }
    for (const auto &entry : remote_table_files) {
        tables.insert(entry.first);
    }

    // Output versions info
    std::ostringstream oss;
    oss << " " << std::left << std::setw(33) << "Table" << std::left << std::setw(7) << "Local"
        << "  Remote\n";
    oss << std::string(1 + 33 + 7 + 2 + 7, '-') << "\n";
    for (const auto &table : tables) {
        int local_version = get_version(local_table_files, table);
        int remote_version = get_version(remote_table_files, table);

        std::string comparator = "==  ";
        if (local_version < remote_version) {
            comparator = "<   ";
        } else if (local_version > remote_version) {
            comparator = ">   ";
        }

        std::string local_version_str =
            local_version == -1 ? "N/A" : "v" + std::to_string(local_version);
        std::string remote_version_str =
            remote_version == -1 ? "N/A" : "v" + std::to_string(remote_version);

        oss << " " << std::left << std::setw(33) << table << std::left << std::setw(5)
            << local_version_str << comparator << remote_version_str << "\n";
    }
    return oss.str();
}

} // namespace pairinteraction
