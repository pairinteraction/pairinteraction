#include "database/Database.hpp"
#include "ket/KetAtom.hpp"
#include "ket/KetAtomCreator.hpp"
#include "utils/paths.hpp"

#include <duckdb.hpp>
#include <fmt/ostream.h>
#include <future>
#include <httplib.h>
#include <nlohmann/json.hpp>
#include <regex>
#include <spdlog/spdlog.h>
#include <sstream>

#if FMT_VERSION < 90000
namespace fmt {
template <typename T>
inline auto streamed(T &&v) {
    return std::forward<T>(v);
}
} // namespace fmt
#endif

Database::Database(bool auto_update)
    : databasedir(paths::get_pairinteraction_cache_directory() / "database"),
      auto_update(auto_update) {

    // Ensure the database directory exists
    if (!std::filesystem::exists(databasedir)) {
        std::filesystem::create_directories(databasedir);
    } else if (!std::filesystem::is_directory(databasedir)) {
        throw std::runtime_error("Database path is not a directory.");
    }

    // Get a dictionary of remotely available databases
    if (auto_update) {
        // Create a pool of clients and call the different endpoints asynchronously
        httplib::Result (httplib::Client::*gf)(const std::string &) = &httplib::Client::Get;
        std::vector<httplib::Client> pool;
        std::vector<std::future<httplib::Result>> futures;

        for (const std::string &endpoint : database_repo_endpoints) {
            pool.emplace_back(httplib::Client("https://api.github.com"));
            pool.back().set_default_headers({{"Accept", "application/vnd.github+json"},
                                             {"X-GitHub-Api-Version", "2022-11-28"}});
            pool.back().set_bearer_token_auth(github_access_token);

            futures.push_back(std::async(std::launch::async, gf, &pool.back(), endpoint));
        }

        // Process the results
        for (auto &future : futures) {
            auto res = future.get();

            // Check if the request was successful
            if (!res) {
                SPDLOG_ERROR("Error accessing database repositories, no response");
                continue;
            } else if (res->status != 200) {
                SPDLOG_ERROR("Error accessing database repositories, status {}", res->status);
                continue;
            }

            // Parse the JSON response
            const std::regex parquet_regex("^(\\w+)_v(\\d+)\\.parquet$");
            auto doc = nlohmann::json::parse(res->body);
            for (auto &asset : doc["assets"]) {
                std::smatch parquet_match;
                auto filename = asset["name"].get<std::string>();
                if (std::regex_match(filename, parquet_match, parquet_regex)) {
                    std::string database = parquet_match[1].str();
                    std::string version_str = parquet_match[2].str();
                    if (std::all_of(version_str.begin(), version_str.end(), ::isdigit)) {
                        int version = std::stoi(version_str);
                        if (remote_databases.count(database) != 0 &&
                            remote_databases.at(database).version >= version) {
                            continue;
                        }
                        remote_databases[database] =
                            RemoteDatabase{asset["url"].get<std::string>().erase(0, 22), version};
                    }
                }
            }
        }
    }
}

template <typename Real>
KetAtom<Real>
Database::get_ket(std::string species, std::optional<Real> energy,
                  std::optional<float> quantum_number_f, std::optional<float> quantum_number_m,
                  std::optional<int> parity, std::optional<int> quantum_number_n,
                  std::optional<Real> quantum_number_nu, std::optional<Real> quantum_number_l,
                  std::optional<Real> quantum_number_s, std::optional<Real> quantum_number_j) {

    ensure_presence_of_local_database(species + "_states");

    // TODO perform database request and delete the following mocked code

    auto ket =
        KetAtom<Real>(energy.value_or(std::numeric_limits<Real>::quiet_NaN()),
                      quantum_number_f.value_or(std::numeric_limits<float>::quiet_NaN()),
                      quantum_number_m.value_or(std::numeric_limits<float>::quiet_NaN()),
                      parity.value_or(-1), "", 1000, species, quantum_number_n.value_or(0),
                      quantum_number_nu.value_or(std::numeric_limits<Real>::quiet_NaN()), 0,
                      quantum_number_l.value_or(std::numeric_limits<Real>::quiet_NaN()), 0,
                      quantum_number_s.value_or(std::numeric_limits<Real>::quiet_NaN()), 0,
                      quantum_number_j.value_or(std::numeric_limits<Real>::quiet_NaN()), 0, *this);

    return ket;
}

template <typename Real>
std::vector<std::shared_ptr<const KetAtom<Real>>> Database::get_kets(
    std::string species, std::optional<Real> min_energy, std::optional<Real> max_energy,
    std::optional<float> min_quantum_number_f, std::optional<float> max_quantum_number_f,
    std::optional<float> min_quantum_number_m, std::optional<float> max_quantum_number_m,
    std::optional<int> parity, std::optional<int> min_quantum_number_n,
    std::optional<int> max_quantum_number_n, std::optional<Real> min_quantum_number_nu,
    std::optional<Real> max_quantum_number_nu, std::optional<Real> min_quantum_number_l,
    std::optional<Real> max_quantum_number_l, std::optional<Real> min_quantum_number_s,
    std::optional<Real> max_quantum_number_s, std::optional<Real> min_quantum_number_j,
    std::optional<Real> max_quantum_number_j) {

    ensure_presence_of_local_database(species + "_states");

    // TODO perform database request and delete the following mocked code

    Real energy = 0;
    float quantum_number_f = 0;
    float quantum_number_m = 0;
    int p = -1;
    size_t id = 1000;
    int quantum_number_n = 0;
    Real quantum_number_nu = 0;
    Real quantum_number_l = 0;
    Real quantum_number_s = 0;
    Real quantum_number_j = 0;

    std::vector<std::shared_ptr<const KetAtom<Real>>> kets;
    kets.push_back(std::make_shared<const KetAtom<Real>>(
        KetAtom<Real>(energy, quantum_number_f, quantum_number_m, p, species, id, species,
                      quantum_number_n, quantum_number_nu, 0, quantum_number_l, 0, quantum_number_s,
                      0, quantum_number_j, 0, *this)));

    return kets;
}

Database::LocalDatabase Database::get_local_database(std::string name) {
    // Loop over all parquet files in the database directory and find the one with the correct name
    // and the highest version
    LocalDatabase db;
    const std::regex parquet_regex("^" + name + "_v(\\d+)\\.parquet$");
    for (const auto &entry : std::filesystem::directory_iterator(databasedir)) {
        if (entry.is_regular_file()) {
            std::smatch parquet_match;
            std::string filename = entry.path().filename().string();
            if (std::regex_match(filename, parquet_match, parquet_regex)) {
                std::string version_str = parquet_match[1].str();
                if (std::all_of(version_str.begin(), version_str.end(), ::isdigit)) {
                    int version = std::stoi(version_str);
                    if (version > db.version) {
                        db.version = version;
                        db.path = entry.path();
                    }
                }
            }
        }
    }
    return db;
}

Database::RemoteDatabase Database::get_remote_database(std::string name) {
    if (remote_databases.count(name) == 0) {
        return RemoteDatabase();
    }
    return remote_databases[name];
}

void Database::ensure_presence_of_local_database(std::string name) {
    // Update the local database if necessary
    if (local_databases.count(name) == 0) {
        auto local_database = get_local_database(name);
        if (auto_update) {
            auto remote_database = get_remote_database(name);
            if (local_database.version < remote_database.version) {
                SPDLOG_INFO("Updating database `{}` from version {} to version {}.", name,
                            local_database.version, remote_database.version);

                httplib::Client cli("https://api.github.com");
                cli.set_default_headers({{"Accept", "application/octet-stream"},
                                         {"X-GitHub-Api-Version", "2022-11-28"}});
                cli.set_bearer_token_auth(github_access_token);
                cli.set_follow_location(true);

                auto res = cli.Get(remote_database.url);
                if (!res || res->status != 200) {
                    SPDLOG_ERROR("Error accessing `{}`: {}", remote_database.url,
                                 fmt::streamed(res.error()));
                } else {
                    if (local_database.version != -1) {
                        std::filesystem::remove(local_database.path);
                    }
                    local_database.version = remote_database.version;
                    local_database.path = databasedir /
                        (name + "_v" + std::to_string(remote_database.version) + ".parquet");
                    std::ofstream out(local_database.path, std::ios::binary);
                    out << res->body;
                    out.close();
                }
            }
        }
        if (local_database.version == -1) {
            if (auto_update) {
                throw std::runtime_error("No database `" + name + "` found.");
            } else {
                throw std::runtime_error("No database `" + name + "` found." +
                                         " Try setting auto_update to true.");
            }
        }
        local_databases[name] = local_database;
    }
}

// Explicit instantiations
template KetAtom<float> Database::get_ket<float>(
    std::string species, std::optional<float> energy, std::optional<float> quantum_number_f,
    std::optional<float> quantum_number_m, std::optional<int> parity,
    std::optional<int> quantum_number_n, std::optional<float> quantum_number_nu,
    std::optional<float> quantum_number_l, std::optional<float> quantum_number_s,
    std::optional<float> quantum_number_j);
template KetAtom<double> Database::get_ket<double>(
    std::string species, std::optional<double> energy, std::optional<float> quantum_number_f,
    std::optional<float> quantum_number_m, std::optional<int> parity,
    std::optional<int> quantum_number_n, std::optional<double> quantum_number_nu,
    std::optional<double> quantum_number_l, std::optional<double> quantum_number_s,
    std::optional<double> quantum_number_j);
template std::vector<std::shared_ptr<const KetAtom<float>>> Database::get_kets<float>(
    std::string species, std::optional<float> min_energy, std::optional<float> max_energy,
    std::optional<float> min_quantum_number_f, std::optional<float> max_quantum_number_f,
    std::optional<float> min_quantum_number_m, std::optional<float> max_quantum_number_m,
    std::optional<int> parity, std::optional<int> min_quantum_number_n,
    std::optional<int> max_quantum_number_n, std::optional<float> min_quantum_number_nu,
    std::optional<float> max_quantum_number_nu, std::optional<float> min_quantum_number_l,
    std::optional<float> max_quantum_number_l, std::optional<float> min_quantum_number_s,
    std::optional<float> max_quantum_number_s, std::optional<float> min_quantum_number_j,
    std::optional<float> max_quantum_number_j);
template std::vector<std::shared_ptr<const KetAtom<double>>> Database::get_kets<double>(
    std::string species, std::optional<double> min_energy, std::optional<double> max_energy,
    std::optional<float> min_quantum_number_f, std::optional<float> max_quantum_number_f,
    std::optional<float> min_quantum_number_m, std::optional<float> max_quantum_number_m,
    std::optional<int> parity, std::optional<int> min_quantum_number_n,
    std::optional<int> max_quantum_number_n, std::optional<double> min_quantum_number_nu,
    std::optional<double> max_quantum_number_nu, std::optional<double> min_quantum_number_l,
    std::optional<double> max_quantum_number_l, std::optional<double> min_quantum_number_s,
    std::optional<double> max_quantum_number_s, std::optional<double> min_quantum_number_j,
    std::optional<double> max_quantum_number_j);
