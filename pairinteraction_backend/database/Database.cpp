#include "database/Database.hpp"
#include "basis/BasisAtom.hpp"
#include "database/AtomDescriptionByParameters.hpp"
#include "database/AtomDescriptionByRanges.hpp"
#include "ket/KetAtom.hpp"
#include "ket/KetAtomCreator.hpp"
#include "operator/OperatorAtom.hpp"
#include "utils/hash.hpp"
#include "utils/ketid.hpp"
#include "utils/paths.hpp"
#include "utils/streamed.hpp"

#include <duckdb.hpp>
#include <fmt/format.h>
#include <fstream>
#include <future>
#include <httplib.h>
#include <nlohmann/json.hpp>
#include <regex>
#include <spdlog/spdlog.h>

Database::Database(bool auto_update)
    : databasedir(paths::get_pairinteraction_cache_directory() / "database"),
      auto_update(auto_update), db(std::make_unique<duckdb::DuckDB>(nullptr)),
      con(std::make_unique<duckdb::Connection>(*db)) {

    const std::regex parquet_regex("^(\\w+)_v(\\d+)\\.parquet$");

    // Ensure the database directory exists
    if (!std::filesystem::exists(databasedir)) {
        std::filesystem::create_directories(databasedir);
    } else if (!std::filesystem::is_directory(databasedir)) {
        throw std::runtime_error("Database path is not a directory.");
    }

    // Ensure that the config directory exists
    std::filesystem::path configdir = paths::get_pairinteraction_config_directory();
    if (!std::filesystem::exists(configdir)) {
        std::filesystem::create_directories(configdir);
    } else if (!std::filesystem::is_directory(configdir)) {
        throw std::runtime_error("Config path is not a directory.");
    }

    // Read in the database_repo_paths if a config file exists, otherwise use the default and
    // write it to the config file
    std::filesystem::path configfile = configdir / "database.json";
    std::string database_repo_host;
    std::vector<std::string> database_repo_paths;
    if (std::filesystem::exists(configfile)) {
        std::ifstream file(configfile);
        nlohmann::json doc;
        file >> doc;

        if (doc.contains("hash") && doc.contains("database_repo_host") &&
            doc.contains("database_repo_paths")) {
            database_repo_host = doc["database_repo_host"].get<std::string>();
            database_repo_paths = doc["database_repo_paths"].get<std::vector<std::string>>();

            // If the values are not equal to the default values but the hash is consistent (i.e.,
            // the user has not changed anything manually), clear the values so that they can be
            // updated
            if (database_repo_host != default_database_repo_host ||
                database_repo_paths != default_database_repo_paths) {
                std::size_t seed = 0;
                hash::hash_combine(seed, database_repo_paths);
                hash::hash_combine(seed, database_repo_host);
                if (seed == doc["hash"].get<std::size_t>()) {
                    database_repo_host.clear();
                    database_repo_paths.clear();
                } else {
                    SPDLOG_INFO("The database repository host and paths have been changed "
                                "manually. Thus, they will not be updated automatically. To reset "
                                "them, delete the file '{}'.",
                                configfile.string());
                }
            }
        }
    }

    // Read in and store the default values if necessary
    if (database_repo_host.empty() || database_repo_paths.empty()) {
        SPDLOG_INFO("Updating the database repository host and paths:");

        database_repo_host = default_database_repo_host;
        database_repo_paths = default_database_repo_paths;
        std::ofstream file(configfile);
        nlohmann::json doc;

        SPDLOG_INFO("* New host: {}", default_database_repo_host);
        SPDLOG_INFO("* New paths: {}", fmt::join(default_database_repo_paths, ", "));

        doc["database_repo_host"] = default_database_repo_host;
        doc["database_repo_paths"] = database_repo_paths;

        std::size_t seed = 0;
        hash::hash_combine(seed, default_database_repo_paths);
        hash::hash_combine(seed, default_database_repo_host);
        doc["hash"] = seed;

        file << doc.dump(4);
    }

    // Create a pool of clients
    if (auto_update) {
        for (size_t i = 0; i < database_repo_paths.size(); i++) {
            pool.emplace_back(httplib::Client(database_repo_host));
            pool.back().set_follow_location(true);
            pool.back().set_connection_timeout(1, 0); // seconds
            pool.back().set_read_timeout(60, 0);      // seconds
            pool.back().set_write_timeout(1, 0);      // seconds
        }
    }

    // Get a dictionary of locally available tables
    for (const auto &entry : std::filesystem::directory_iterator(databasedir)) {
        if (entry.is_regular_file()) {
            std::smatch parquet_match;
            std::string filename = entry.path().filename().string();
            if (std::regex_match(filename, parquet_match, parquet_regex)) {
                std::string name = parquet_match[1].str();
                std::string version_str = parquet_match[2].str();
                if (std::all_of(version_str.begin(), version_str.end(), ::isdigit)) {
                    int version = std::stoi(version_str);
                    if (version > tables[name].local_version) {
                        tables[name].local_version = version;
                        tables[name].local_path = entry.path();
                    }
                }
            }
        }
    }

    // Get a dictionary of remotely available tables
    if (auto_update) {
        // Call the different endpoints asynchronously
        httplib::Result (httplib::Client::*gf)(const std::string &, const httplib::Headers &) =
            &httplib::Client::Get;
        std::vector<std::future<httplib::Result>> futures;
        std::vector<std::filesystem::path> filenames;

        for (size_t i = 0; i < database_repo_paths.size(); i++) {
            // Get the last modified date of the last JSON response
            std::string lastmodified = "";
            filenames.push_back(databasedir /
                                ("latest_" +
                                 std::to_string(std::hash<std::string>{}(database_repo_paths[i])) +
                                 ".json"));
            if (std::filesystem::exists(filenames.back())) {
                std::ifstream file(filenames.back());
                nlohmann::json doc;
                file >> doc;
                if (doc.contains("last-modified")) {
                    lastmodified = doc["last-modified"].get<std::string>();
                }
            }

            // Create a get request conditioned on the last-modified date
            httplib::Headers headers;
            if (lastmodified.empty()) {
                headers = {{"X-GitHub-Api-Version", "2022-11-28"},
                           {"Accept", "application/vnd.github+json"}};
            } else {
                headers = {{"X-GitHub-Api-Version", "2022-11-28"},
                           {"Accept", "application/vnd.github+json"},
                           {"if-modified-since", lastmodified}};
            }
            futures.push_back(
                std::async(std::launch::async, gf, &pool[i], database_repo_paths[i], headers));
        }

        // Process the results
        for (size_t i = 0; i < database_repo_paths.size(); i++) {
            SPDLOG_INFO("Accessing database repository path: {}", database_repo_paths[i]);
            auto res = futures[i].get();

            // Check if the request was successful
            if (!res) {
                if (res.error() != httplib::Error::Unknown ||
                    !std::filesystem::exists(filenames[i])) {
                    SPDLOG_ERROR("Access error: {}", fmt::streamed(res.error()));
                    continue;
                }

            } else if ((res->status == 403 || res->status == 429) &&
                       res->has_header("x-ratelimit-remaining") &&
                       res->get_header_value("x-ratelimit-remaining") == "0") {
                int waittime = -1;
                if (res->has_header("retry-after")) {
                    waittime = std::stoi(res->get_header_value("retry-after"));
                } else if (res->has_header("x-ratelimit-reset")) {
                    waittime =
                        std::stoi(res->get_header_value("x-ratelimit-reset")) - time(nullptr);
                }
                SPDLOG_ERROR("Access error, rate limit exceeded. Auto update "
                             "is disabled. You should not retry until after {} seconds.",
                             waittime);
                auto_update = false;
                continue;

            } else if (res->status != 200) {
                SPDLOG_ERROR("Access error, status {}.", res->status);
                continue;
            }

            // Parse and store the JSON response together with the last-modified date
            nlohmann::json doc;
            if (!res || res->status == 304) {
                SPDLOG_INFO("Using cached overview of available tables.");
                std::ifstream file(filenames[i]);
                file >> doc;
            } else {
                SPDLOG_INFO("Using downloaded overview of available tables.");
                doc = nlohmann::json::parse(res->body);
                if (res->has_header("last-modified")) {
                    doc["last-modified"] = res->get_header_value("last-modified");
                }
                std::ofstream file(filenames[i]);
                file << doc;
            }

            // Process the assets
            for (auto &asset : doc["assets"]) {
                std::smatch parquet_match;
                auto filename = asset["name"].get<std::string>();
                if (std::regex_match(filename, parquet_match, parquet_regex)) {
                    std::string name = parquet_match[1].str();
                    std::string version_str = parquet_match[2].str();
                    if (std::all_of(version_str.begin(), version_str.end(), ::isdigit)) {
                        int version = std::stoi(version_str);
                        if (version > tables[name].remote_version) {
                            tables[name].remote_version = version;
                            tables[name].remote_path = asset["url"].get<std::string>().erase(0, 22);
                        }
                    }
                }
            }
        }
    }

    // Print availability of tables
    auto species_availability = get_availability_of_species();
    auto wigner_availability = get_availability_of_wigner_table();
    SPDLOG_INFO("Availability of database tables for species and Wigner 3j symbols:");
    for (const auto &s : species_availability) {
        SPDLOG_INFO("* {} (locally available: {}, up to date: {}, fully downloaded: {})", s.name,
                    s.locally_available, s.up_to_date, s.fully_downloaded);
    }
    SPDLOG_INFO("* Wigner 3j symbols (locally available: {}, up to date: {})",
                wigner_availability.locally_available, wigner_availability.up_to_date);
}

Database::~Database() = default;

std::vector<Database::AvailabilitySpecies> Database::get_availability_of_species() {
    std::vector<AvailabilitySpecies> availability;

    // Get a list of all available species
    for (const auto &[name, table] : tables) {
        // Ensure that the table is a states table
        if (name.size() < 7 || name.substr(name.size() - 7) != "_states") {
            continue;
        }

        // Read off the species name
        std::string species_name = name.substr(0, name.size() - 7);

        // Append the species
        availability.push_back({species_name, table.local_version != -1,
                                table.local_version >= table.remote_version,
                                table.local_version != -1});
    }

    // Check whether all tables are downloaded and the downloaded tables up to date
    std::vector<std::string> identifier_of_tables = {"states",
                                                     "matrix_elements_d",
                                                     "matrix_elements_q",
                                                     "matrix_elements_o",
                                                     "matrix_elements_mu",
                                                     "matrix_elements_dia"};
    for (size_t i = 0; i < availability.size(); i++) {
        for (const auto &identifier : identifier_of_tables) {
            std::string name = availability[i].name + "_" + identifier;
            if (tables.count(name) > 0) {
                if (tables[name].local_version == -1) {
                    availability[i].fully_downloaded = false;
                } else if (tables[name].local_version < tables[name].remote_version) {
                    availability[i].up_to_date = false;
                }
            }
        }
    }

    return availability;
}

Database::AvailabilityWigner Database::get_availability_of_wigner_table() {
    std::string name = "wigner";
    if (tables.count(name) == 0) {
        return {false, false};
    }
    return {tables[name].local_version != -1,
            tables[name].local_version >= tables[name].remote_version};
}

template <typename Real>
KetAtom<Real> Database::get_ket(std::string species,
                                const AtomDescriptionByParameters<Real> &description) {

    ensure_presence_of_table(species + "_states");

    // Check that the specifications are valid
    if (description.quantum_number_n.has_value()) {
        ensure_quantum_number_n_is_allowed(species + "_states");
    }
    if (!description.quantum_number_m.has_value() &&
        description.quantum_number_f.value_or(1) != 0) {
        throw std::runtime_error("The quantum number m must be specified if f is not zero.");
    }
    if (description.quantum_number_f.has_value() &&
        2 * description.quantum_number_f.value() !=
            std::rintf(2 * description.quantum_number_f.value())) {
        throw std::runtime_error("The quantum number f must be an integer or half-integer.");
    }
    if (description.quantum_number_f.has_value() && description.quantum_number_f.value() < 0) {
        throw std::runtime_error("The quantum number f must be positive.");
    }
    if (description.quantum_number_m.has_value() &&
        2 * description.quantum_number_m.value() !=
            std::rintf(2 * description.quantum_number_m.value())) {
        throw std::runtime_error("The quantum number m must be an integer or half-integer.");
    }

    // Describe the state
    std::string where = "";
    std::string separator = "";
    if (description.energy.has_value()) {
        // The following condition derives from demanding that quantum number n that corresponds to
        // the energy "E_n = -1/(2*n^2)" is not off by more than 1 from the actual quantum number n,
        // i.e., "sqrt(-1/(2*E_n)) - sqrt(-1/(2*E_{n-1})) = 1"
        where += separator +
            fmt::format("SQRT(-1/(2*energy)) BETWEEN {} AND {}",
                        std::sqrt(-1 / (2 * description.energy.value())) - 0.5,
                        std::sqrt(-1 / (2 * description.energy.value())) + 0.5);
        separator = " AND ";
    }
    if (description.quantum_number_f.has_value()) {
        where += separator + fmt::format("f = {}", description.quantum_number_f.value());
        separator = " AND ";
    }
    if (description.parity.has_value()) {
        where += separator + fmt::format("parity = {}", description.parity.value());
        separator = " AND ";
    }
    if (description.quantum_number_n.has_value()) {
        where += separator + fmt::format("n = {}", description.quantum_number_n.value());
        separator = " AND ";
    }
    if (description.quantum_number_nu.has_value()) {
        where += separator +
            fmt::format("exp_nu BETWEEN {} AND {}", description.quantum_number_nu.value() - 0.5,
                        description.quantum_number_nu.value() + 0.5);
        separator = " AND ";
    }
    if (description.quantum_number_l.has_value()) {
        where += separator +
            fmt::format("exp_l BETWEEN {} AND {}", description.quantum_number_l.value() - 0.5,
                        description.quantum_number_l.value() + 0.5);
        separator = " AND ";
    }
    if (description.quantum_number_s.has_value()) {
        where += separator +
            fmt::format("exp_s BETWEEN {} AND {}", description.quantum_number_s.value() - 0.5,
                        description.quantum_number_s.value() + 0.5);
        separator = " AND ";
    }
    if (description.quantum_number_j.has_value()) {
        where += separator +
            fmt::format("exp_j BETWEEN {} AND {}", description.quantum_number_j.value() - 0.5,
                        description.quantum_number_j.value() + 0.5);
        separator = " AND ";
    }
    if (separator.empty()) {
        where += "FALSE";
    }

    std::string orderby = "";
    separator = "";
    if (description.energy.has_value()) {
        orderby += separator +
            fmt::format("(SQRT(-1/(2*energy)) - {})^2",
                        std::sqrt(-1 / (2 * description.energy.value())));
        separator = " + ";
    }
    if (description.quantum_number_nu.has_value()) {
        orderby +=
            separator + fmt::format("(exp_nu - {})^2", description.quantum_number_nu.value());
        separator = " + ";
    }
    if (description.quantum_number_l.has_value()) {
        orderby += separator + fmt::format("(exp_l - {})^2", description.quantum_number_l.value());
        separator = " + ";
    }
    if (description.quantum_number_s.has_value()) {
        orderby += separator + fmt::format("(exp_s - {})^2", description.quantum_number_s.value());
        separator = " + ";
    }
    if (description.quantum_number_j.has_value()) {
        orderby += separator + fmt::format("(exp_j - {})^2", description.quantum_number_j.value());
        separator = " + ";
    }
    if (separator.empty()) {
        orderby += "id";
    }

    // Ask the database for the described state
    auto result = con->Query(fmt::format(
        R"(SELECT energy, f, parity, id, n, exp_nu, std_nu, exp_l, std_l, exp_s, std_s,
        exp_j, std_j FROM '{}' WHERE {} ORDER BY {} ASC LIMIT 1)",
        tables[species + "_states"].local_path.string(), where, orderby));

    if (result->HasError()) {
        throw std::runtime_error("Error querying the database: " + result->GetError());
    }

    if (result->RowCount() == 0) {
        throw std::runtime_error("No state found.");
    }

    // Check the types of the columns
    const auto &types = result->types;
    const auto &labels = result->names;
    const std::vector<duckdb::LogicalType> ref_types = {
        duckdb::LogicalType::DOUBLE, duckdb::LogicalType::DOUBLE, duckdb::LogicalType::BIGINT,
        duckdb::LogicalType::BIGINT, duckdb::LogicalType::BIGINT, duckdb::LogicalType::DOUBLE,
        duckdb::LogicalType::DOUBLE, duckdb::LogicalType::DOUBLE, duckdb::LogicalType::DOUBLE,
        duckdb::LogicalType::DOUBLE, duckdb::LogicalType::DOUBLE, duckdb::LogicalType::DOUBLE,
        duckdb::LogicalType::DOUBLE};

    for (size_t i = 0; i < types.size(); i++) {
        if (types[i] != ref_types[i]) {
            throw std::runtime_error("Wrong type for '" + labels[i] + "'.");
        }
    }
    // Construct the state
    auto chunk = result->Fetch();

    auto result_quantum_number_m = description.quantum_number_m.value_or(0);
    auto result_energy = duckdb::FlatVector::GetData<double>(chunk->data[0])[0];
    auto result_quantum_number_f = duckdb::FlatVector::GetData<double>(chunk->data[1])[0];
    auto result_parity = duckdb::FlatVector::GetData<int64_t>(chunk->data[2])[0];
    auto result_id = ketid::atom::get(duckdb::FlatVector::GetData<int64_t>(chunk->data[3])[0],
                                      result_quantum_number_m);
    auto result_quantum_number_n = duckdb::FlatVector::GetData<int64_t>(chunk->data[4])[0];
    auto result_quantum_number_nu_exp = duckdb::FlatVector::GetData<double>(chunk->data[5])[0];
    auto result_quantum_number_nu_std = duckdb::FlatVector::GetData<double>(chunk->data[6])[0];
    auto result_quantum_number_l_exp = duckdb::FlatVector::GetData<double>(chunk->data[7])[0];
    auto result_quantum_number_l_std = duckdb::FlatVector::GetData<double>(chunk->data[8])[0];
    auto result_quantum_number_s_exp = duckdb::FlatVector::GetData<double>(chunk->data[9])[0];
    auto result_quantum_number_s_std = duckdb::FlatVector::GetData<double>(chunk->data[10])[0];
    auto result_quantum_number_j_exp = duckdb::FlatVector::GetData<double>(chunk->data[11])[0];
    auto result_quantum_number_j_std = duckdb::FlatVector::GetData<double>(chunk->data[12])[0];

    return KetAtom<Real>(result_energy, result_quantum_number_f, result_quantum_number_m,
                         result_parity, result_id, species, result_quantum_number_n,
                         result_quantum_number_nu_exp, result_quantum_number_nu_std,
                         result_quantum_number_l_exp, result_quantum_number_l_std,
                         result_quantum_number_s_exp, result_quantum_number_s_std,
                         result_quantum_number_j_exp, result_quantum_number_j_std);
}

template <typename Scalar>
BasisAtom<Scalar> Database::get_basis(std::string species,
                                      const AtomDescriptionByRanges<Scalar> &description,
                                      std::vector<size_t> additional_ket_ids) {
    using real_t = typename traits::NumTraits<Scalar>::real_t;

    ensure_presence_of_table(species + "_states");

    // Check that the specifications are valid
    if (description.min_quantum_number_n.has_value() ||
        description.max_quantum_number_n.has_value()) {
        ensure_quantum_number_n_is_allowed(species + "_states");
    }

    // Describe the states
    std::string where = "(";
    std::string separator = "";
    if (description.min_energy.has_value()) {
        where += separator + fmt::format("energy >= {}", description.min_energy.value());
        separator = " AND ";
    }
    if (description.max_energy.has_value()) {
        where += separator + fmt::format("energy <= {}", description.max_energy.value());
        separator = " AND ";
    }
    if (description.min_quantum_number_f.has_value()) {
        where += separator + fmt::format("f >= {}", description.min_quantum_number_f.value());
        separator = " AND ";
    }
    if (description.max_quantum_number_f.has_value()) {
        where += separator + fmt::format("f <= {}", description.max_quantum_number_f.value());
        separator = " AND ";
    }
    if (description.min_quantum_number_m.has_value()) {
        where += separator + fmt::format("m >= {}", description.min_quantum_number_m.value());
        separator = " AND ";
    }
    if (description.max_quantum_number_m.has_value()) {
        where += separator + fmt::format("m <= {}", description.max_quantum_number_m.value());
        separator = " AND ";
    }
    if (description.parity.has_value()) {
        where += separator + fmt::format("parity = {}", description.parity.value());
        separator = " AND ";
    }
    if (description.min_quantum_number_n.has_value()) {
        where += separator + fmt::format("n >= {}", description.min_quantum_number_n.value());
        separator = " AND ";
    }
    if (description.max_quantum_number_n.has_value()) {
        where += separator + fmt::format("n <= {}", description.max_quantum_number_n.value());
        separator = " AND ";
    }
    if (description.min_quantum_number_nu.has_value()) {
        where += separator +
            fmt::format("exp_nu >= {}-2*std_nu", description.min_quantum_number_nu.value());
        separator = " AND ";
    }
    if (description.max_quantum_number_nu.has_value()) {
        where += separator +
            fmt::format("exp_nu <= {}+2*std_nu", description.max_quantum_number_nu.value());
        separator = " AND ";
    }
    if (description.min_quantum_number_l.has_value()) {
        where += separator +
            fmt::format("exp_l >= {}-2*std_l", description.min_quantum_number_l.value());
        separator = " AND ";
    }
    if (description.max_quantum_number_l.has_value()) {
        where += separator +
            fmt::format("exp_l <= {}+2*std_l", description.max_quantum_number_l.value());
        separator = " AND ";
    }
    if (description.min_quantum_number_s.has_value()) {
        where += separator +
            fmt::format("exp_s >= {}-2*std_s", description.min_quantum_number_s.value());
        separator = " AND ";
    }
    if (description.max_quantum_number_s.has_value()) {
        where += separator +
            fmt::format("exp_s <= {}+2*std_s", description.max_quantum_number_s.value());
        separator = " AND ";
    }
    if (description.min_quantum_number_j.has_value()) {
        where += separator +
            fmt::format("exp_j >= {}-2*std_j", description.min_quantum_number_j.value());
        separator = " AND ";
    }
    if (description.max_quantum_number_j.has_value()) {
        where += separator +
            fmt::format("exp_j <= {}+2*std_j", description.max_quantum_number_j.value());
        separator = " AND ";
    }
    if (separator.empty()) {
        where += "FALSE";
    }
    where += ")";
    if (!additional_ket_ids.empty()) {
        where += fmt::format(" OR {} IN ({})", ketid::atom::SQL_TERM,
                             fmt::join(additional_ket_ids, ","));
    }

    // Create a table containing the described states
    std::string uuid;
    {
        auto result = con->Query(R"(SELECT UUID()::varchar)");
        if (result->HasError()) {
            throw std::runtime_error("Error selecting uuid: " + result->GetError());
        }
        uuid =
            duckdb::FlatVector::GetData<duckdb::string_t>(result->Fetch()->data[0])[0].GetString();
    }
    {
        auto result = con->Query(fmt::format(
            R"(CREATE TABLE '{}' AS SELECT *, {} AS ketid FROM (
                SELECT *,
                UNNEST(list_transform(generate_series(0,(2*f)::bigint),
                x -> x::double-f)) AS m FROM '{}'
            ) WHERE {})",
            uuid, ketid::atom::SQL_TERM, tables[species + "_states"].local_path.string(), where));

        if (result->HasError()) {
            throw std::runtime_error("Error creating table: " + result->GetError());
        }
    }

    // Ask the table for the extreme values of the quantum numbers
    {
        std::string select = "";
        std::string separator = "";
        if (description.min_energy.has_value()) {
            select += separator + fmt::format("MIN(energy) AS min_energy");
            separator = ", ";
        }
        if (description.max_energy.has_value()) {
            select += separator + fmt::format("MAX(energy) AS max_energy");
            separator = ", ";
        }
        if (description.min_quantum_number_f.has_value()) {
            select += separator + fmt::format("MIN(f) AS min_f");
            separator = ", ";
        }
        if (description.max_quantum_number_f.has_value()) {
            select += separator + fmt::format("MAX(f) AS max_f");
            separator = ", ";
        }
        if (description.min_quantum_number_m.has_value()) {
            select += separator + fmt::format("MIN(m) AS min_m");
            separator = ", ";
        }
        if (description.max_quantum_number_m.has_value()) {
            select += separator + fmt::format("MAX(m) AS max_m");
            separator = ", ";
        }
        if (description.min_quantum_number_n.has_value()) {
            select += separator + fmt::format("MIN(n) AS min_n");
            separator = ", ";
        }
        if (description.max_quantum_number_n.has_value()) {
            select += separator + fmt::format("MAX(n) AS max_n");
            separator = ", ";
        }
        if (description.min_quantum_number_nu.has_value()) {
            select += separator + fmt::format("MIN(exp_nu) AS min_nu");
            separator = ", ";
        }
        if (description.max_quantum_number_nu.has_value()) {
            select += separator + fmt::format("MAX(exp_nu) AS max_nu");
            separator = ", ";
        }
        if (description.min_quantum_number_l.has_value()) {
            select += separator + fmt::format("MIN(exp_l) AS min_l");
            separator = ", ";
        }
        if (description.max_quantum_number_l.has_value()) {
            select += separator + fmt::format("MAX(exp_l) AS max_l");
            separator = ", ";
        }
        if (description.min_quantum_number_s.has_value()) {
            select += separator + fmt::format("MIN(exp_s) AS min_s");
            separator = ", ";
        }
        if (description.max_quantum_number_s.has_value()) {
            select += separator + fmt::format("MAX(exp_s) AS max_s");
            separator = ", ";
        }
        if (description.min_quantum_number_j.has_value()) {
            select += separator + fmt::format("MIN(exp_j) AS min_j");
            separator = ", ";
        }
        if (description.max_quantum_number_j.has_value()) {
            select += separator + fmt::format("MAX(exp_j) AS max_j");
            separator = ", ";
        }

        if (!separator.empty()) {
            auto result = con->Query(fmt::format(R"(SELECT {} FROM '{}')", select, uuid));

            if (result->HasError()) {
                throw std::runtime_error("Error querying the database: " + result->GetError());
            }

            auto chunk = result->Fetch();

            for (size_t i = 0; i < chunk->ColumnCount(); i++) {
                if (duckdb::FlatVector::IsNull(chunk->data[i], 0)) {
                    throw std::runtime_error("No state found.");
                }
            }

            size_t idx = 0;
            if (description.min_energy.has_value()) {
                auto min_energy = duckdb::FlatVector::GetData<double>(chunk->data[idx++])[0];
                if (std::sqrt(-1 / (2 * min_energy)) - 1 >
                    std::sqrt(-1 / (2 * description.min_energy.value()))) {
                    SPDLOG_WARN("No state found with the requested minimum energy. Requested: {}, "
                                "found: {}.",
                                description.min_energy.value(), min_energy);
                }
            }
            if (description.max_energy.has_value()) {
                auto max_energy = duckdb::FlatVector::GetData<double>(chunk->data[idx++])[0];
                if (std::sqrt(-1 / (2 * max_energy)) + 1 <
                    std::sqrt(-1 / (2 * description.max_energy.value()))) {
                    SPDLOG_WARN("No state found with the requested maximum energy. Requested: {}, "
                                "found: {}.",
                                description.max_energy.value(), max_energy);
                }
            }
            if (description.min_quantum_number_f.has_value()) {
                auto min_f = duckdb::FlatVector::GetData<double>(chunk->data[idx++])[0];
                if (min_f > description.min_quantum_number_f.value()) {
                    SPDLOG_WARN("No state found with the requested minimum quantum number f. "
                                "Requested: {}, found: {}.",
                                description.min_quantum_number_f.value(), min_f);
                }
            }
            if (description.max_quantum_number_f.has_value()) {
                auto max_f = duckdb::FlatVector::GetData<double>(chunk->data[idx++])[0];
                if (max_f < description.max_quantum_number_f.value()) {
                    SPDLOG_WARN("No state found with the requested maximum quantum number f. "
                                "Requested: {}, found: {}.",
                                description.max_quantum_number_f.value(), max_f);
                }
            }
            if (description.min_quantum_number_m.has_value()) {
                auto min_m = duckdb::FlatVector::GetData<double>(chunk->data[idx++])[0];
                if (min_m > description.min_quantum_number_m.value()) {
                    SPDLOG_WARN("No state found with the requested minimum quantum number m. "
                                "Requested: {}, found: {}.",
                                description.min_quantum_number_m.value(), min_m);
                }
            }
            if (description.max_quantum_number_m.has_value()) {
                auto max_m = duckdb::FlatVector::GetData<double>(chunk->data[idx++])[0];
                if (max_m < description.max_quantum_number_m.value()) {
                    SPDLOG_WARN("No state found with the requested maximum quantum number m. "
                                "Requested: {}, found: {}.",
                                description.max_quantum_number_m.value(), max_m);
                }
            }
            if (description.min_quantum_number_n.has_value()) {
                auto min_n = duckdb::FlatVector::GetData<int64_t>(chunk->data[idx++])[0];
                if (min_n > description.min_quantum_number_n.value()) {
                    SPDLOG_WARN("No state found with the requested minimum quantum number n. "
                                "Requested: {}, found: {}.",
                                description.min_quantum_number_n.value(), min_n);
                }
            }
            if (description.max_quantum_number_n.has_value()) {
                auto max_n = duckdb::FlatVector::GetData<int64_t>(chunk->data[idx++])[0];
                if (max_n < description.max_quantum_number_n.value()) {
                    SPDLOG_WARN("No state found with the requested maximum quantum number n. "
                                "Requested: {}, found: {}.",
                                description.max_quantum_number_n.value(), max_n);
                }
            }
            if (description.min_quantum_number_nu.has_value()) {
                auto min_nu = duckdb::FlatVector::GetData<double>(chunk->data[idx++])[0];
                if (min_nu - 1 > description.min_quantum_number_nu.value()) {
                    SPDLOG_WARN("No state found with the requested minimum quantum number nu. "
                                "Requested: {}, found: {}.",
                                description.min_quantum_number_nu.value(), min_nu);
                }
            }
            if (description.max_quantum_number_nu.has_value()) {
                auto max_nu = duckdb::FlatVector::GetData<double>(chunk->data[idx++])[0];
                if (max_nu + 1 < description.max_quantum_number_nu.value()) {
                    SPDLOG_WARN("No state found with the requested maximum quantum number nu. "
                                "Requested: {}, found: {}.",
                                description.max_quantum_number_nu.value(), max_nu);
                }
            }
            if (description.min_quantum_number_l.has_value()) {
                auto min_l = duckdb::FlatVector::GetData<double>(chunk->data[idx++])[0];
                if (min_l - 1 > description.min_quantum_number_l.value()) {
                    SPDLOG_WARN("No state found with the requested minimum quantum number l. "
                                "Requested: {}, found: {}.",
                                description.min_quantum_number_l.value(), min_l);
                }
            }
            if (description.max_quantum_number_l.has_value()) {
                auto max_l = duckdb::FlatVector::GetData<double>(chunk->data[idx++])[0];
                if (max_l + 1 < description.max_quantum_number_l.value()) {
                    SPDLOG_WARN("No state found with the requested maximum quantum number l. "
                                "Requested: {}, found: {}.",
                                description.max_quantum_number_l.value(), max_l);
                }
            }
            if (description.min_quantum_number_s.has_value()) {
                auto min_s = duckdb::FlatVector::GetData<double>(chunk->data[idx++])[0];
                if (min_s - 1 > description.min_quantum_number_s.value()) {
                    SPDLOG_WARN("No state found with the requested minimum quantum number s. "
                                "Requested: {}, found: {}.",
                                description.min_quantum_number_s.value(), min_s);
                }
            }
            if (description.max_quantum_number_s.has_value()) {
                auto max_s = duckdb::FlatVector::GetData<double>(chunk->data[idx++])[0];
                if (max_s + 1 < description.max_quantum_number_s.value()) {
                    SPDLOG_WARN("No state found with the requested maximum quantum number s. "
                                "Requested: {}, found: {}.",
                                description.max_quantum_number_s.value(), max_s);
                }
            }
            if (description.min_quantum_number_j.has_value()) {
                auto min_j = duckdb::FlatVector::GetData<double>(chunk->data[idx++])[0];
                if (min_j - 1 > description.min_quantum_number_j.value()) {
                    SPDLOG_WARN("No state found with the requested minimum quantum number j. "
                                "Requested: {}, found: {}.",
                                description.min_quantum_number_j.value(), min_j);
                }
            }
            if (description.max_quantum_number_j.has_value()) {
                auto max_j = duckdb::FlatVector::GetData<double>(chunk->data[idx++])[0];
                if (max_j + 1 < description.max_quantum_number_j.value()) {
                    SPDLOG_WARN("No state found with the requested maximum quantum number j. "
                                "Requested: {}, found: {}.",
                                description.max_quantum_number_j.value(), max_j);
                }
            }
        }
    }

    // Ask the table for the described states
    auto result = con->Query(fmt::format(
        R"(SELECT energy, f, m, parity, ketid, n, exp_nu, std_nu, exp_l, std_l,
        exp_s, std_s, exp_j, std_j FROM '{}' ORDER BY ketid ASC)",
        uuid));

    if (result->HasError()) {
        throw std::runtime_error("Error querying the database: " + result->GetError());
    }

    if (result->RowCount() == 0) {
        throw std::runtime_error("No state found.");
    }

    // Check the types of the columns
    const auto &types = result->types;
    const auto &labels = result->names;
    const std::vector<duckdb::LogicalType> ref_types = {
        duckdb::LogicalType::DOUBLE, duckdb::LogicalType::DOUBLE, duckdb::LogicalType::DOUBLE,
        duckdb::LogicalType::BIGINT, duckdb::LogicalType::BIGINT, duckdb::LogicalType::BIGINT,
        duckdb::LogicalType::DOUBLE, duckdb::LogicalType::DOUBLE, duckdb::LogicalType::DOUBLE,
        duckdb::LogicalType::DOUBLE, duckdb::LogicalType::DOUBLE, duckdb::LogicalType::DOUBLE,
        duckdb::LogicalType::DOUBLE, duckdb::LogicalType::DOUBLE};

    for (size_t i = 0; i < types.size(); i++) {
        if (types[i] != ref_types[i]) {
            throw std::runtime_error("Wrong type for '" + labels[i] + "'.");
        }
    }

    // Construct the states
    std::vector<std::shared_ptr<const KetAtom<real_t>>> kets;
    kets.reserve(result->RowCount());
    double last_energy = std::numeric_limits<double>::lowest();

    for (auto chunk = result->Fetch(); chunk; chunk = result->Fetch()) {

        auto chunk_energy = duckdb::FlatVector::GetData<double>(chunk->data[0]);
        auto chunk_quantum_number_f = duckdb::FlatVector::GetData<double>(chunk->data[1]);
        auto chunk_quantum_number_m = duckdb::FlatVector::GetData<double>(chunk->data[2]);
        auto chunk_parity = duckdb::FlatVector::GetData<int64_t>(chunk->data[3]);
        auto chunk_id = duckdb::FlatVector::GetData<int64_t>(chunk->data[4]);
        auto chunk_quantum_number_n = duckdb::FlatVector::GetData<int64_t>(chunk->data[5]);
        auto chunk_quantum_number_nu_exp = duckdb::FlatVector::GetData<double>(chunk->data[6]);
        auto chunk_quantum_number_nu_std = duckdb::FlatVector::GetData<double>(chunk->data[7]);
        auto chunk_quantum_number_l_exp = duckdb::FlatVector::GetData<double>(chunk->data[8]);
        auto chunk_quantum_number_l_std = duckdb::FlatVector::GetData<double>(chunk->data[9]);
        auto chunk_quantum_number_s_exp = duckdb::FlatVector::GetData<double>(chunk->data[10]);
        auto chunk_quantum_number_s_std = duckdb::FlatVector::GetData<double>(chunk->data[11]);
        auto chunk_quantum_number_j_exp = duckdb::FlatVector::GetData<double>(chunk->data[12]);
        auto chunk_quantum_number_j_std = duckdb::FlatVector::GetData<double>(chunk->data[13]);

        for (size_t i = 0; i < chunk->size(); i++) {

            // Check that the states are sorted by energy
            if (chunk_energy[i] < last_energy) {
                throw std::runtime_error("The states are not sorted by energy.");
            }
            last_energy = chunk_energy[i];

            // Append a new state
            kets.push_back(std::make_shared<const KetAtom<real_t>>(KetAtom<real_t>(
                chunk_energy[i], chunk_quantum_number_f[i], chunk_quantum_number_m[i],
                chunk_parity[i], chunk_id[i], species, chunk_quantum_number_n[i],
                chunk_quantum_number_nu_exp[i], chunk_quantum_number_nu_std[i],
                chunk_quantum_number_l_exp[i], chunk_quantum_number_l_std[i],
                chunk_quantum_number_s_exp[i], chunk_quantum_number_s_std[i],
                chunk_quantum_number_j_exp[i], chunk_quantum_number_j_std[i])));
        }
    }

    return BasisAtom<Scalar>(std::move(kets), uuid, *this, species);
}

template <typename Scalar>
OperatorAtom<Scalar> Database::get_operator(const BasisAtom<Scalar> &basis, OperatorType type,
                                            int q) {
    std::string specifier;
    int kappa;
    switch (type) {
    case OperatorType::DIPOLE:
        specifier = "matrix_elements_d";
        kappa = 1;
        break;
    case OperatorType::QUADRUPOLE:
        specifier = "matrix_elements_q";
        kappa = 2;
        break;
    case OperatorType::OCTUPOLE:
        specifier = "matrix_elements_o";
        kappa = 3;
        break;
    case OperatorType::MAGNETICDIPOLE:
        specifier = "matrix_elements_mu";
        kappa = 1;
        break;
    case OperatorType::DIAMAGNETIC:
        specifier = "matrix_elements_dia";
        kappa = 0;
        break;
    case OperatorType::ENERGY:
        specifier = "energy";
        kappa = 0;
        break;
    default:
        throw std::runtime_error("Unknown operator type.");
    }

    ensure_presence_of_table("wigner");
    if (specifier != "energy") {
        ensure_presence_of_table(basis.species + "_" + specifier);
    }

    // Check that the specifications are valid
    if (std::abs(q) > kappa) {
        throw std::runtime_error("Invalid q.");
    }

    // Ask the database for the operator
    duckdb::unique_ptr<duckdb::MaterializedQueryResult> result;
    if (specifier != "energy") {
        result = con->Query(fmt::format(
            R"(WITH s AS (
                SELECT id, f, m, ketid FROM '{}'
            ),
            b AS (
                SELECT MIN(f) AS min_f, MAX(f) AS max_f,
                MIN(id) AS min_id, MAX(id) AS max_id
                FROM s
            ),
            w_filtered AS (
                SELECT *
                FROM '{}'
                WHERE kappa = {} AND q = {} AND
                f_initial BETWEEN (SELECT min_f FROM b) AND (SELECT max_f FROM b) AND
                f_final BETWEEN (SELECT min_f FROM b) AND (SELECT max_f FROM b)
            ),
            e_filtered AS (
                SELECT *
                FROM '{}'
                WHERE
                id_initial BETWEEN (SELECT min_id FROM b) AND (SELECT max_id FROM b) AND
                id_final BETWEEN (SELECT min_id FROM b) AND (SELECT max_id FROM b)
            )
            SELECT
            s2.ketid AS row,
            s1.ketid AS col,
            e.val*w.val AS val
            FROM e_filtered AS e
            JOIN s AS s1 ON e.id_initial = s1.id
            JOIN s AS s2 ON e.id_final = s2.id
            JOIN w_filtered AS w ON
            w.f_initial = s1.f AND w.m_initial = s1.m AND
            w.f_final = s2.f AND w.m_final = s2.m
            ORDER BY row ASC, col ASC)",
            basis.table, tables["wigner"].local_path.string(), kappa, q,
            tables[basis.species + "_" + specifier].local_path.string()));
    } else {
        result = con->Query(fmt::format(
            R"(SELECT ketid as row, ketid as col, energy as val FROM '{}' ORDER BY row ASC)",
            basis.table));
    }

    if (result->HasError()) {
        throw std::runtime_error("Error querying the database: " + result->GetError());
    }

    if (result->RowCount() == 0) {
        throw std::runtime_error("No matrix elements founds.");
    }

    // Check the types of the columns
    const auto &types = result->types;
    const auto &labels = result->names;
    const std::vector<duckdb::LogicalType> ref_types = {
        duckdb::LogicalType::BIGINT, duckdb::LogicalType::BIGINT, duckdb::LogicalType::DOUBLE};
    for (size_t i = 0; i < types.size(); i++) {
        if (types[i] != ref_types[i]) {
            throw std::runtime_error("Wrong type for '" + labels[i] + "'.");
        }
    }

    // Construct the matrix
    int dim = basis.get_number_of_states();
    int num_entries = result->RowCount();

    std::vector<int> outerIndexPtr;
    std::vector<int> innerIndices;
    std::vector<Scalar> values;
    outerIndexPtr.reserve(dim + 1);
    innerIndices.reserve(num_entries);
    values.reserve(num_entries);

    int last_row = -1;

    for (auto chunk = result->Fetch(); chunk; chunk = result->Fetch()) {

        auto chunk_row = duckdb::FlatVector::GetData<int64_t>(chunk->data[0]);
        auto chunk_col = duckdb::FlatVector::GetData<int64_t>(chunk->data[1]);
        auto chunk_val = duckdb::FlatVector::GetData<double>(chunk->data[2]);

        for (size_t i = 0; i < chunk->size(); i++) {
            int row = basis.ket_id_to_index.at(chunk_row[i]);
            if (row != last_row) {
                if (row < last_row) {
                    throw std::runtime_error("The rows are not sorted.");
                }
                for (; last_row < row; last_row++) {
                    outerIndexPtr.push_back(innerIndices.size());
                }
            }
            innerIndices.push_back(basis.ket_id_to_index.at(chunk_col[i]));
            values.push_back(chunk_val[i]);
        }
    }

    for (; last_row < dim + 1; last_row++) {
        outerIndexPtr.push_back(innerIndices.size());
    }

    Eigen::Map<const Eigen::SparseMatrix<Scalar, Eigen::RowMajor>> matrix_map(
        dim, dim, num_entries, outerIndexPtr.data(), innerIndices.data(), values.data());

    // Transform the matrix into the provided basis and return it
    Eigen::SparseMatrix<Scalar, Eigen::RowMajor> matrix =
        basis.coefficients.adjoint() * matrix_map * basis.coefficients;

    // Construct the operator and return it
    return OperatorAtom(basis, type, q, std::move(matrix));
}

void Database::ensure_presence_of_table(std::string name) {
    if (tables.count(name) == 0 && auto_update) {
        throw std::runtime_error("No database `" + name + "` found.");
    }

    if (tables.count(name) == 0 && !auto_update) {
        throw std::runtime_error("No database `" + name +
                                 "` found. Try setting auto_update to true.");
    }

    if (auto_update && tables[name].local_version < tables[name].remote_version) {
        SPDLOG_INFO("Updating database `{}` from version {} to version {}.", name,
                    tables[name].local_version, tables[name].remote_version);
        auto res = pool.front().Get(
            tables[name].remote_path.string(),
            {{"X-GitHub-Api-Version", "2022-11-28"}, {"Accept", "application/octet-stream"}});
        if (!res || res->status != 200) {
            SPDLOG_ERROR("Error accessing `{}`: {}", tables[name].remote_path.string(),
                         fmt::streamed(res.error()));
        } else if ((res->status == 403 || res->status == 429) &&
                   res->has_header("x-ratelimit-remaining") &&
                   res->get_header_value("x-ratelimit-remaining") == "0") {
            int waittime = -1;
            if (res->has_header("retry-after")) {
                waittime = std::stoi(res->get_header_value("retry-after"));
            } else if (res->has_header("x-ratelimit-reset")) {
                waittime = std::stoi(res->get_header_value("x-ratelimit-reset")) - time(nullptr);
            }
            SPDLOG_ERROR("Error accessing database repositories, rate limit exceeded. Auto update "
                         "is disabled. You should not retry until after {} seconds.",
                         waittime);
            auto_update = false;
        } else {
            if (tables[name].local_version != -1) {
                std::filesystem::remove(tables[name].local_path);
            }
            tables[name].local_version = tables[name].remote_version;
            tables[name].local_path = databasedir /
                (name + "_v" + std::to_string(tables[name].remote_version) + ".parquet");
            std::ofstream out(tables[name].local_path, std::ios::binary);
            out << res->body;
            out.close();
        }
    }
}

void Database::ensure_quantum_number_n_is_allowed(std::string name) {
    auto result =
        con->Query(fmt::format(R"(SELECT n FROM '{}' LIMIT 1)", tables[name].local_path.string()));

    if (result->HasError()) {
        throw std::runtime_error("Error querying the database: " + result->GetError());
    }

    if (result->RowCount() == 0) {
        throw std::runtime_error("No state found.");
    }

    auto chunk = result->Fetch();

    if (chunk->data[0].GetType() != duckdb::LogicalType::BIGINT) {
        throw std::runtime_error("Wrong type for 'n'.");
    }

    if (duckdb::FlatVector::GetData<int64_t>(chunk->data[0])[0] <= 0) {
        throw std::runtime_error(
            "The specified species does not have a well-defined principal quantum number n. "
            "Use the effective principal quantum number nu instead.");
    }
}

Database &Database::get_global_instance() {
    thread_local static Database database;
    return database;
}

// Explicit instantiations
template KetAtom<float>
Database::get_ket<float>(std::string species,
                         const AtomDescriptionByParameters<float> &description);
template KetAtom<double>
Database::get_ket<double>(std::string species,
                          const AtomDescriptionByParameters<double> &description);
template BasisAtom<float>
Database::get_basis<float>(std::string species, const AtomDescriptionByRanges<float> &description,
                           std::vector<size_t> additional_ket_ids);
template BasisAtom<double>
Database::get_basis<double>(std::string species, const AtomDescriptionByRanges<double> &description,
                            std::vector<size_t> additional_ket_ids);
template BasisAtom<std::complex<float>> Database::get_basis<std::complex<float>>(
    std::string species, const AtomDescriptionByRanges<std::complex<float>> &description,
    std::vector<size_t> additional_ket_ids);
template BasisAtom<std::complex<double>> Database::get_basis<std::complex<double>>(
    std::string species, const AtomDescriptionByRanges<std::complex<double>> &description,
    std::vector<size_t> additional_ket_ids);
template OperatorAtom<float> Database::get_operator<float>(const BasisAtom<float> &basis,
                                                           OperatorType type, int q);
template OperatorAtom<double> Database::get_operator<double>(const BasisAtom<double> &basis,
                                                             OperatorType type, int q);
template OperatorAtom<std::complex<float>>
Database::get_operator<std::complex<float>>(const BasisAtom<std::complex<float>> &basis,
                                            OperatorType type, int q);
template OperatorAtom<std::complex<double>>
Database::get_operator<std::complex<double>>(const BasisAtom<std::complex<double>> &basis,
                                             OperatorType type, int q);

///////////////////////////////////////////////////////////////////////////////////////
// Test cases
///////////////////////////////////////////////////////////////////////////////////////

#include "utils/streamed.hpp"
#include <doctest/doctest.h>
#include <spdlog/spdlog.h>

DOCTEST_TEST_CASE("get a KetAtom") {
    Database &database = Database::get_global_instance();

    AtomDescriptionByParameters<float> description;
    description.quantum_number_n = 60;
    description.quantum_number_l = 0;
    description.quantum_number_m = 0;

    auto ket = database.get_ket<float>("Rb", description);

    SPDLOG_LOGGER_INFO(spdlog::get("doctest"), "KetAtom: {}", fmt::streamed(ket));
}

DOCTEST_TEST_CASE("get a BasisAtom") {
    Database &database = Database::get_global_instance();

    AtomDescriptionByRanges<float> description;
    description.min_quantum_number_n = 60;
    description.max_quantum_number_n = 60;
    description.min_quantum_number_l = 0;
    description.max_quantum_number_l = 1;

    auto basis = database.get_basis<float>("Rb", description, {});

    for (const auto &ket : basis) {
        SPDLOG_LOGGER_INFO(spdlog::get("doctest"), "KetAtom: {}", fmt::streamed(ket));
    }
}

DOCTEST_TEST_CASE("get an OperatorAtom") {
    Database &database = Database::get_global_instance();

    AtomDescriptionByRanges<float> description;
    description.min_quantum_number_n = 60;
    description.max_quantum_number_n = 60;
    description.min_quantum_number_l = 0;
    description.max_quantum_number_l = 1;

    auto basis = database.get_basis<float>("Rb", description, {});

    auto dipole = database.get_operator<float>(basis, OperatorType::DIPOLE, 0);

    SPDLOG_LOGGER_INFO(spdlog::get("doctest"), "Number of basis states: {}",
                       basis.get_number_of_states());
    SPDLOG_LOGGER_INFO(spdlog::get("doctest"), "Number of non-zero entries: {}",
                       dipole.get_matrix().nonZeros());
}
