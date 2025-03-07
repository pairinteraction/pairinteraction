#include "pairinteraction/database/ParquetManager.hpp"

#include "pairinteraction/database/GitHubDownloader.hpp"
#include "pairinteraction/version.hpp"

#include <cpptrace/cpptrace.hpp>
#include <ctime>
#include <duckdb.hpp>
#include <filesystem>
#include <fmt/core.h>
#include <fstream>
#include <future>
#include <iomanip>
#include <miniz.h>
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
    if (!std::filesystem::exists(directory_ / "tables")) {
        std::filesystem::create_directories(directory_ / "tables");
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
    remote_asset_info.clear();

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

        // Update remote_asset_info based on the asset entries
        for (auto &asset : doc["assets"]) {
            std::string name = asset["name"].get<std::string>();
            std::smatch match;

            if (std::regex_match(name, match, remote_regex) && match.size() == 4) {
                std::string key = match[1].str();
                int version_major = std::stoi(match[2].str());
                int version_minor = std::stoi(match[3].str());

                if (version_major != COMPATIBLE_DATABASE_VERSION_MAJOR) {
                    throw std::runtime_error("Pairinteraction is incompatible with the database "
                                             "repositories. Consider upgrading pairinteraction.");
                }

                auto it = remote_asset_info.find(key);
                if (it == remote_asset_info.end() || version_minor > it->second.version_minor) {
                    std::string remote_url = asset["url"].get<std::string>();
                    const std::string host = downloader.get_host();
                    remote_asset_info[key] = {version_minor, remote_url.erase(0, host.size())};
                }
            }
        }
    }
}

void ParquetManager::scan_local() {
    local_asset_info.clear();

    // Iterate over files in the directory to update local_asset_info
    for (const auto &entry : std::filesystem::directory_iterator(directory_ / "tables")) {
        std::string name = entry.path().filename().string();
        std::smatch match;

        if (entry.is_directory() && std::regex_match(name, match, local_regex) &&
            match.size() == 4) {
            std::string key = match[1].str();
            int version_major = std::stoi(match[2].str());
            int version_minor = std::stoi(match[3].str());

            if (version_major != COMPATIBLE_DATABASE_VERSION_MAJOR) {
                continue;
            }

            auto it = local_asset_info.find(key);
            if (it == local_asset_info.end() || version_minor > it->second.version_minor) {
                local_asset_info[key].version_minor = version_minor;
                for (const auto &subentry : std::filesystem::directory_iterator(entry)) {
                    if (subentry.is_regular_file() && subentry.path().extension() == ".parquet") {
                        local_asset_info[key].paths[subentry.path().stem().string()] = {
                            subentry.path().string(), false};
                    }
                }
            }
        }
    }
}

void ParquetManager::react_on_rate_limit_reached(std::time_t reset_time) {
    repo_paths_.clear();
    remote_asset_info.clear();
    SPDLOG_WARN("Rate limit reached, resets at {}. The download of database tables is disabled.",
                format_time(reset_time));
}

void ParquetManager::update_local_asset(const std::string &key) {
    // Get remote version if available
    int remote_version = -1;
    auto remote_it = remote_asset_info.find(key);
    if (remote_it != remote_asset_info.end()) {
        remote_version = remote_it->second.version_minor;
    }

    // Get local version if available and check if it is up-to-date
    {
        int local_version = -1;
        std::shared_lock<std::shared_mutex> lock(mtx_local);
        auto local_it = local_asset_info.find(key);
        if (local_it != local_asset_info.end()) {
            local_version = local_it->second.version_minor;
        }
        if (local_version >= remote_version) {
            return;
        }
    }

    // If it is not up-to-date, acquire a unique lock for updating the table
    std::unique_lock<std::shared_mutex> lock(mtx_local);

    // Re-check if the table is up to date because another thread might have updated it
    int local_version = -1;
    auto local_it = local_asset_info.find(key);
    if (local_it != local_asset_info.end()) {
        local_version = local_it->second.version_minor;
    }
    if (local_version >= remote_version) {
        return;
    }

    // Download the remote file
    std::string endpoint = remote_it->second.endpoint;
    SPDLOG_INFO("Downloading {}_v{}.{} from {}", key, COMPATIBLE_DATABASE_VERSION_MAJOR,
                remote_version, endpoint);

    auto result = downloader.download(endpoint, "", true).get();
    if (result.status_code == 403 || result.status_code == 429) {
        react_on_rate_limit_reached(result.rate_limit.reset_time);
        return;
    }
    if (result.status_code != 200) {
        throw std::runtime_error(fmt::format("Failed to download table {}: status code {}.",
                                             endpoint, result.status_code));
    }

    // Unzip the downloaded file
    mz_zip_archive zip_archive{};
    mz_bool status =
        mz_zip_reader_init_mem(&zip_archive, result.body.data(), result.body.size(), 0);
    if (status == 0) {
        throw std::runtime_error("Failed to initialize zip archive.");
    }

    for (mz_uint i = 0; i < mz_zip_reader_get_num_files(&zip_archive); i++) {
        mz_zip_archive_file_stat file_stat;
        if (mz_zip_reader_file_stat(&zip_archive, i, &file_stat) == 0) {
            throw std::runtime_error("Failed to get file stat from zip archive.");
        }

        // Skip directories
        const char *filename = static_cast<const char *>(file_stat.m_filename);
        size_t len = std::strlen(filename);
        if (len > 0 && filename[len - 1] == '/') {
            continue;
        }

        auto path = directory_ / "tables" / file_stat.m_filename;
        SPDLOG_INFO("Storing table to {}", path.string());

        // Extract the file to memory
        std::vector<char> buffer(file_stat.m_uncomp_size);
        if (mz_zip_reader_extract_to_mem(&zip_archive, i, buffer.data(), buffer.size(), 0) == 0) {
            throw std::runtime_error(fmt::format("Failed to extract {}.", file_stat.m_filename));
        }

        // Ensure the parent directory exists
        std::filesystem::create_directories(path.parent_path());

        // Save the extracted file
        std::ofstream out(path.string(), std::ios::binary);
        if (!out) {
            throw std::runtime_error(fmt::format("Failed to open {} for writing", path.string()));
        }
        out.write(buffer.data(), buffer.size());
        out.close();

        // Update the local asset/table info
        local_asset_info[key].version_minor = remote_version;
        local_asset_info[key].paths[path.stem().string()] = {path.string(), false};
    }

    mz_zip_reader_end(&zip_archive);
}

void ParquetManager::cache_table(std::unordered_map<std::string, PathInfo>::iterator table_it) {
    // Check if the table is already cached
    {
        std::shared_lock<std::shared_mutex> lock(mtx_local);
        if (table_it->second.cached) {
            return;
        }
    }

    // Acquire a unique lock for caching the table
    std::unique_lock<std::shared_mutex> lock(mtx_local);

    // Re-check if the table is already cached because another thread might have cached it
    if (table_it->second.cached) {
        return;
    }

    // Cache the table in memory
    std::string table_name;
    {
        auto result = con.Query(R"(SELECT UUID()::varchar)");
        if (result->HasError()) {
            throw cpptrace::runtime_error("Error selecting a unique table name: " +
                                          result->GetError());
        }
        table_name =
            duckdb::FlatVector::GetData<duckdb::string_t>(result->Fetch()->data[0])[0].GetString();
    }

    {
        auto result = con.Query(fmt::format(R"(CREATE TEMP TABLE '{}' AS SELECT * FROM '{}')",
                                            table_name, table_it->second.path));
        if (result->HasError()) {
            throw cpptrace::runtime_error("Error creating table: " + result->GetError());
        }
    }

    table_it->second.path = table_name;
    table_it->second.cached = true;
}

std::string ParquetManager::get_path(const std::string &key, const std::string &table) {
    // Update the local table if a newer version is available remotely
    this->update_local_asset(key);

    // Ensure availability of the local table file
    auto asset_it = local_asset_info.find(key);
    if (asset_it == local_asset_info.end()) {
        throw std::runtime_error("Table " + key + "_" + table + " not found.");
    }
    auto table_it = asset_it->second.paths.find(table);
    if (table_it == asset_it->second.paths.end()) {
        throw std::runtime_error("Table " + key + "_" + table + " not found.");
    }

    // Cache the local table in memory if requested
    if (use_cache_) {
        this->cache_table(table_it);
    }

    // Return the path to the local table file
    return table_it->second.path;
}

std::string ParquetManager::get_versions_info() const {
    // Helper lambda returns the version string if available
    auto get_version = [](const auto &map, const std::string &table) -> int {
        if (auto it = map.find(table); it != map.end()) {
            return it->second.version_minor;
        }
        return -1;
    };

    // Gather all unique table names
    std::set<std::string> tables;
    for (const auto &entry : local_asset_info) {
        tables.insert(entry.first);
    }
    for (const auto &entry : remote_asset_info) {
        tables.insert(entry.first);
    }

    // Output versions info
    std::ostringstream oss;
    oss << " " << std::left << std::setw(33) << "Asset" << std::left << std::setw(8) << "Local"
        << "  Remote\n";
    oss << std::string(1 + 33 + 8 + 2 + 7, '-') << "\n";
    for (const auto &table : tables) {
        int local_version = get_version(local_asset_info, table);
        int remote_version = get_version(remote_asset_info, table);

        std::string comparator = "==  ";
        if (local_version < remote_version) {
            comparator = "<   ";
        } else if (local_version > remote_version) {
            comparator = ">   ";
        }

        std::string local_version_str = local_version == -1
            ? "N/A"
            : "v" + std::to_string(COMPATIBLE_DATABASE_VERSION_MAJOR) + "." +
                std::to_string(local_version);
        std::string remote_version_str = remote_version == -1
            ? "N/A"
            : "v" + std::to_string(COMPATIBLE_DATABASE_VERSION_MAJOR) + "." +
                std::to_string(remote_version);

        oss << " " << std::left << std::setw(33) << table << std::left << std::setw(6)
            << local_version_str << comparator << remote_version_str << "\n";
    }
    return oss.str();
}

} // namespace pairinteraction
