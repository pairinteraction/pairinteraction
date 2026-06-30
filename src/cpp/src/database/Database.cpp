// SPDX-FileCopyrightText: 2024 PairInteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "pairinteraction/database/Database.hpp"

#include "pairinteraction/basis/BasisAtom.hpp"
#include "pairinteraction/database/AtomDescriptionByParameters.hpp"
#include "pairinteraction/database/AtomDescriptionByRanges.hpp"
#include "pairinteraction/database/GitHubDownloader.hpp"
#include "pairinteraction/database/ParquetManager.hpp"
#include "pairinteraction/enums/OperatorType.hpp"
#include "pairinteraction/ket/KetAtom.hpp"
#include "pairinteraction/utils/TaskControl.hpp"
#include "pairinteraction/utils/hash.hpp"
#include "pairinteraction/utils/ket_id.hpp"
#include "pairinteraction/utils/paths.hpp"
#include "pairinteraction/utils/streamed.hpp"

#include <algorithm>
#include <cpptrace/cpptrace.hpp>
#include <duckdb.hpp>
#include <fmt/core.h>
#include <fmt/format.h>
#include <fmt/ranges.h>
#include <fstream>
#include <iterator>
#include <nlohmann/json.hpp>
#include <oneapi/tbb.h>
#include <spdlog/spdlog.h>
#include <string>
#include <system_error>
#include <unordered_map>
#include <unordered_set>

namespace pairinteraction {

namespace {

std::string format_expectation_value_range(const std::string &value_column,
                                           const std::string &std_column,
                                           const Range<double> &range,
                                           double standard_deviation_factor) {
    return fmt::format("{} BETWEEN {}-{}*{} AND {}+{}*{}", value_column, range.min(),
                       standard_deviation_factor, std_column, range.max(),
                       standard_deviation_factor, std_column);
}

// Find the index of the result column with the given name.
size_t get_column_index(const std::vector<std::string> &names, const std::string &name) {
    auto it = std::find(names.begin(), names.end(), name);
    if (it == names.end()) {
        throw std::runtime_error("Missing database column '" + name + "'.");
    }
    return static_cast<size_t>(std::distance(names.begin(), it));
}

// Read a single value of a duckdb result column as a double, regardless of its logical type.
double get_entry_as_double(duckdb::Vector &vector, const duckdb::LogicalType &type, size_t row) {
    switch (type.id()) {
    case duckdb::LogicalTypeId::DOUBLE:
        return duckdb::FlatVector::GetData<double>(vector)[row];
    case duckdb::LogicalTypeId::BIGINT:
        return static_cast<double>(duckdb::FlatVector::GetData<int64_t>(vector)[row]);
    case duckdb::LogicalTypeId::BOOLEAN:
        return duckdb::FlatVector::GetData<bool>(vector)[row] ? 1.0 : 0.0;
    default:
        throw std::runtime_error("Cannot read database column of type " + type.ToString() +
                                 " as a quantum number.");
    }
}

struct QuantumNumbers {
    std::unordered_map<std::string, double> values;
    std::unordered_map<std::string, double> stds;
};

QuantumNumbers get_quantum_numbers_from_row(duckdb::DataChunk &chunk,
                                            const std::vector<duckdb::LogicalType> &types,
                                            const std::vector<std::string> &names,
                                            const std::unordered_set<std::string> &excluded_columns,
                                            size_t row) {
    QuantumNumbers quantum_numbers;
    for (size_t col = 0; col < names.size(); ++col) {
        const std::string &name = names[col];
        if (excluded_columns.contains(name)) {
            continue;
        }
        double value = get_entry_as_double(chunk.data[col], types[col], row);
        if (name.starts_with("std_")) {
            quantum_numbers.stds[name.substr(4)] = value;
        } else if (name.starts_with("exp_")) {
            quantum_numbers.values[name.substr(4)] = value;
        } else {
            quantum_numbers.values[name] = value;
        }
    }
    return quantum_numbers;
}
} // namespace

void ensure_consistent_quantum_numbers(double quantum_number_f, double quantum_number_m) {
    if (2 * quantum_number_m != std::rint(2 * quantum_number_m)) {
        throw std::runtime_error("The quantum number m must be an integer or half-integer.");
    }
    if (2 * quantum_number_f != std::rint(2 * quantum_number_f)) {
        throw std::runtime_error("The quantum number f must be an integer or half-integer.");
    }
    if (quantum_number_f + quantum_number_m != std::rint(quantum_number_f + quantum_number_m)) {
        throw std::invalid_argument(
            "The quantum numbers f and m must be both either integers or half-integers.");
    }
    if (std::abs(quantum_number_m) > quantum_number_f) {
        throw std::invalid_argument(
            "The absolute value of the quantum number m must be less than or equal to f.");
    }
}

Database::Database() : Database(default_download_missing) {}

Database::Database(bool download_missing)
    : Database(download_missing, default_use_cache, default_database_dir) {}

Database::Database(std::filesystem::path database_dir)
    : Database(default_download_missing, default_use_cache, std::move(database_dir)) {}

Database::Database(bool download_missing, bool use_cache, std::filesystem::path database_dir)
    : download_missing_(download_missing), use_cache_(use_cache),
      database_dir_(std::move(database_dir)), db(std::make_unique<duckdb::DuckDB>(nullptr)),
      con(std::make_unique<duckdb::Connection>(*db)) {

    if (database_dir_.empty()) {
        database_dir_ = default_database_dir;
    }

    // Ensure the database directory exists
    if (!std::filesystem::exists(database_dir_)) {
        std::filesystem::create_directories(database_dir_);
    }
    database_dir_ = std::filesystem::canonical(database_dir_);
    if (!std::filesystem::is_directory(database_dir_)) {
        throw std::filesystem::filesystem_error("Cannot access database", database_dir_.string(),
                                                std::make_error_code(std::errc::not_a_directory));
    }
    SPDLOG_INFO("Using database directory: {}", database_dir_.string());

    // Ensure that the config directory exists
    std::filesystem::path configdir = paths::get_config_directory();
    if (!std::filesystem::exists(configdir)) {
        std::filesystem::create_directories(configdir);
    } else if (!std::filesystem::is_directory(configdir)) {
        throw std::filesystem::filesystem_error("Cannot access config directory ",
                                                configdir.string(),
                                                std::make_error_code(std::errc::not_a_directory));
    }

    // Read in the database_repo_paths if a config file exists, otherwise use the default and
    // write it to the config file
    std::filesystem::path configfile = configdir / "database.json";
    std::string database_repo_host;
    std::vector<std::string> database_repo_paths;
    if (std::filesystem::exists(configfile)) {
        std::ifstream file(configfile);
        nlohmann::json doc = nlohmann::json::parse(file, nullptr, false);

        if (!doc.is_discarded() && doc.contains("hash") && doc.contains("database_repo_host") &&
            doc.contains("database_repo_paths")) {
            database_repo_host = doc["database_repo_host"].get<std::string>();
            database_repo_paths = doc["database_repo_paths"].get<std::vector<std::string>>();

            // If the values are not equal to the default values but the hash is consistent (i.e.,
            // the user has not changed anything manually), clear the values so that they can be
            // updated
            if (database_repo_host != default_database_repo_host ||
                database_repo_paths != default_database_repo_paths) {
                std::size_t seed = 0;
                utils::hash_combine(seed, database_repo_paths);
                utils::hash_combine(seed, database_repo_host);
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
        utils::hash_combine(seed, default_database_repo_paths);
        utils::hash_combine(seed, default_database_repo_host);
        doc["hash"] = seed;

        file << doc.dump(4);
    }

    // Limit the memory usage of duckdb's buffer manager
    {
        auto result = con->Query("PRAGMA max_memory = '8GB';");
        if (result->HasError()) {
            throw cpptrace::runtime_error("Error setting the memory limit: " + result->GetError());
        }
    }

    // Instantiate a database manager that provides access to database tables. If a table
    // is outdated/not available locally, it will be downloaded if download_missing_ is true.
    if (!download_missing_) {
        database_repo_paths.clear();
    }
    downloader = std::make_unique<GitHubDownloader>();
    manager = std::make_unique<ParquetManager>(database_dir_, *downloader, database_repo_paths,
                                               *con, use_cache_);
    manager->scan_local();
    manager->scan_remote();

    // Print versions of tables
    std::istringstream iss(manager->get_versions_info());
    for (std::string line; std::getline(iss, line);) {
        SPDLOG_INFO(line);
    }
}

Database::~Database() = default;

const std::unordered_set<std::string> &Database::get_column_names(const std::string &table_path) {
    if (auto it = column_names_cache.find(table_path); it != column_names_cache.end()) {
        return it->second;
    }
    auto result = con->Query(fmt::format(R"(SELECT * FROM '{}' LIMIT 0)", table_path));
    if (result->HasError()) {
        throw cpptrace::runtime_error("Error querying the database columns: " + result->GetError());
    }
    std::unordered_set<std::string> names(result->names.begin(), result->names.end());

    // Every states table must provide the columns required for constructing kets and must not
    // contain an 'm' column, since m added separately below.
    for (const auto &required : {"id", "f", "energy"}) {
        if (!names.contains(required)) {
            throw std::runtime_error(
                fmt::format("The database table '{}' is missing the required column '{}'.",
                            table_path, required));
        }
    }
    if (names.contains("m")) {
        throw std::runtime_error(
            fmt::format("The database table '{}' must not contain a column 'm'.", table_path));
    }

    // If another thread inserted the same entry concurrently, insert keeps the existing value and
    // returns an iterator to it, so the reference stays valid (elements are never erased).
    return column_names_cache.insert({table_path, std::move(names)}).first->second;
}

std::shared_ptr<const KetAtom> Database::get_ket(const std::string &species,
                                                 const AtomDescriptionByParameters &description) {
    // Check that the specifications are valid
    if (!description.quantum_numbers.contains("m")) {
        throw std::invalid_argument("The quantum number m must be specified.");
    }
    for (const auto &[name, value] : description.quantum_numbers) {
        if ((name == "f" || name == "m") && 2 * value != std::rint(2 * value)) {
            throw std::invalid_argument("The quantum number " + name +
                                        " must be an integer or half-integer.");
        }
        if (name == "f" && value < 0) {
            throw std::invalid_argument("The quantum number " + name + " must be positive.");
        }
    }

    const auto &columns = get_column_names(manager->get_path(species, "states"));

    // Describe the state. The quantum numbers n, f and parity are matched exactly, while all other
    // quantum numbers are matched within a +-0.5 window (they can deviate from the requested value,
    // e.g. expectation values in MQDT). The result is ordered by the distance to the requested
    // values, so the nearest state is returned.
    std::string where;
    std::string where_separator;
    std::string orderby;
    std::string orderby_separator;
    if (description.energy.has_value()) {
        // The following condition derives from demanding that quantum number n that corresponds to
        // the energy "E_n = -1/(2*n^2)" is not off by more than 1 from the actual quantum number n,
        // i.e., "sqrt(-1/(2*E_n)) - sqrt(-1/(2*E_{n-1})) = 1"
        double n_from_energy = std::sqrt(-1 / (2 * description.energy.value()));
        where += where_separator +
            fmt::format("SQRT(-1/(2*energy)) BETWEEN {} AND {}", n_from_energy - 0.5,
                        n_from_energy + 0.5);
        where_separator = " AND ";
        orderby += orderby_separator + fmt::format("(SQRT(-1/(2*energy)) - {})^2", n_from_energy);
        orderby_separator = " + ";
    }
    for (const auto &[name, value] : description.quantum_numbers) {
        if (name == "m") {
            continue; // m is encoded into the id, not stored as a queryable column
        }
        std::string column = columns.contains("exp_" + name) ? "exp_" + name : name;
        if (!columns.contains(column)) {
            throw std::invalid_argument(
                fmt::format("The quantum number '{}' is not stored in the database table for "
                            "species '{}'.",
                            name, species));
        }
        double tolerance = (name == "n" || name == "f" || name == "parity") ? 0.0 : 0.5;
        where += where_separator +
            fmt::format("{} BETWEEN {} AND {}", column, value - tolerance, value + tolerance);
        where_separator = " AND ";
        orderby += orderby_separator + fmt::format("({} - {})^2", column, value);
        orderby_separator = " + ";
    }
    if (where_separator.empty()) {
        where += "FALSE";
    }
    if (orderby_separator.empty()) {
        orderby += "id";
    }

    // Ask the database for the described state
    set_task_status("Loading atomic ket from database...");
    auto result = con->Query(fmt::format(
        R"(SELECT *, {} AS order_val FROM '{}' WHERE {} ORDER BY order_val ASC LIMIT 2)", orderby,
        manager->get_path(species, "states"), where));

    if (result->HasError()) {
        throw cpptrace::runtime_error("Error querying the database: " + result->GetError());
    }

    if (result->RowCount() == 0) {
        throw std::invalid_argument("No state found.");
    }

    // Get the first chunk of the results (the first chunk is sufficient as we need two rows at
    // most). Every column except energy, id and the synthetic order_val is treated as a quantum
    // number; m is not a database column and is injected from the description.
    const auto &types = result->types;
    const auto &names = result->names;
    const std::unordered_set<std::string> excluded_columns = {"energy", "id", "order_val"};
    auto chunk = result->Fetch();

    size_t energy_column = get_column_index(names, "energy");
    size_t id_column = get_column_index(names, "id");
    double quantum_number_m = description.quantum_numbers.at("m");

    auto make_ket = [&](size_t row) {
        auto quantum_numbers =
            get_quantum_numbers_from_row(*chunk, types, names, excluded_columns, row);
        quantum_numbers.values["m"] = quantum_number_m;
        double energy = get_entry_as_double(chunk->data[energy_column], types[energy_column], row);
        auto id =
            utils::encode_as_ket_id({.id = static_cast<size_t>(duckdb::FlatVector::GetData<int64_t>(
                                         chunk->data[id_column])[row]),
                                     .m = quantum_number_m});
        return KetAtom(typename KetAtom::Private(), energy, species,
                       std::move(quantum_numbers.values), std::move(quantum_numbers.stds), *this,
                       id);
    };

    // Describe a candidate row by its species and quantum numbers (used for error messages)
    auto describe_ket = [&](size_t row) {
        auto quantum_numbers =
            get_quantum_numbers_from_row(*chunk, types, names, excluded_columns, row);
        std::string description = species;
        for (const auto &[name, value] : quantum_numbers.values) {
            description += fmt::format(" {}={}", name, value);
        }
        description += fmt::format(" m={}", quantum_number_m);
        return description;
    };

    // Check that the ket is uniquely specified
    if (chunk->size() > 1) {
        size_t order_val_column = get_column_index(names, "order_val");
        auto order_val_0 =
            get_entry_as_double(chunk->data[order_val_column], types[order_val_column], 0);
        auto order_val_1 =
            get_entry_as_double(chunk->data[order_val_column], types[order_val_column], 1);

        if (order_val_1 - order_val_0 <= order_val_0) {
            throw std::invalid_argument(
                fmt::format("The ket is not uniquely specified. Possible kets are:\n{}\n{}",
                            describe_ket(0), describe_ket(1)));
        }
    }

    // Construct the state
    auto ket = make_ket(0);

    // Check database consistency
    ensure_consistent_quantum_numbers(ket.get_quantum_number("f"), quantum_number_m);

    return std::make_shared<const KetAtom>(std::move(ket));
}

template <typename Scalar>
std::shared_ptr<const BasisAtom<Scalar>>
Database::get_basis(const std::string &species, const AtomDescriptionByRanges &description,
                    const std::vector<size_t> &additional_ket_ids) {
    // The quantum number m is restricted separately because it is generated by UNNEST below.
    auto range_quantum_number_m = [&description]() {
        auto it = description.quantum_number_ranges.find("m");
        return it != description.quantum_number_ranges.end() ? it->second : Range<double>{};
    }();

    const auto &columns = get_column_names(manager->get_path(species, "states"));

    // Describe the states by all restrictions that do not involve the quantum number m
    std::string where = "(";
    std::string separator;
    if (description.range_energy.is_finite()) {
        where += separator +
            fmt::format("energy BETWEEN {} AND {}", description.range_energy.min(),
                        description.range_energy.max());
        separator = " AND ";
    }
    for (const auto &[name, range] : description.quantum_number_ranges) {
        if (name == "m" || !range.is_finite()) {
            continue;
        }
        std::string exp_column = "exp_" + name;
        std::string std_column = "std_" + name;
        if (columns.contains(exp_column) && columns.contains(std_column)) {
            where += separator +
                format_expectation_value_range(
                         exp_column, std_column, range,
                         description.quantum_number_standard_deviation_factor);
        } else if (columns.contains(name)) {
            where +=
                separator + fmt::format("{} BETWEEN {} AND {}", name, range.min(), range.max());
        } else {
            throw std::invalid_argument(
                fmt::format("The quantum number '{}' is not stored in the database table for "
                            "species '{}'.",
                            name, species));
        }
        separator = " AND ";
    }
    if (separator.empty()) {
        // If the description contains no restrictions at all, it describes no states
        where += range_quantum_number_m.is_finite() ? "TRUE" : "FALSE";
    }
    where += ")";

    // Describe the restriction of the quantum number m
    std::string where_m = "(";
    if (range_quantum_number_m.is_finite()) {
        where_m += fmt::format("m BETWEEN {} AND {}", range_quantum_number_m.min(),
                               range_quantum_number_m.max());
    } else {
        where_m += "TRUE";
    }
    where_m += ")";

    // Create a table containing the described states
    std::string canonical_basis_id;
    {
        auto result = con->Query(R"(SELECT UUID()::varchar)");
        if (result->HasError()) {
            throw cpptrace::runtime_error("Error selecting canonical_basis_id: " +
                                          result->GetError());
        }
        canonical_basis_id =
            duckdb::FlatVector::GetData<duckdb::string_t>(result->Fetch()->data[0])[0].GetString();
    }
    {
        set_task_status("Selecting atomic basis states...");
        auto result = con->Query(fmt::format(
            R"(CREATE TEMP TABLE '{}' AS SELECT *, id*{}+(2*m+{})::bigint AS ketid FROM (
                SELECT *,
                UNNEST(list_transform(generate_series(0,(2*f)::bigint),
                x -> x::double-f)) AS m FROM (
                    SELECT * FROM '{}' WHERE {}
                )
            ) WHERE {})",
            canonical_basis_id, utils::KET_ID_STRIDE, utils::M_OFFSET,
            manager->get_path(species, "states"), where, where_m));

        if (result->HasError()) {
            throw cpptrace::runtime_error("Error creating table: " + result->GetError());
        }
    }
    {
        auto result = con->Query(
            fmt::format(R"(ALTER TABLE '{}' ADD PRIMARY KEY (ketid))", canonical_basis_id));

        if (result->HasError()) {
            throw cpptrace::runtime_error("Error adding primary key: " + result->GetError());
        }
    }

    // Add the additional kets to the table if they are not already contained in it
    if (!additional_ket_ids.empty()) {
        set_task_status("Selecting additional kets...");
        std::vector<std::string> values;
        values.reserve(additional_ket_ids.size());
        for (size_t ket_id : additional_ket_ids) {
            auto [id, m] = utils::decode_from_ket_id(ket_id);
            values.push_back(fmt::format("({}, {}::double, {})", id, m, ket_id));
        }
        auto result = con->Query(fmt::format(
            R"(INSERT OR IGNORE INTO '{0}'
            SELECT s.*, v.m, v.ketid FROM '{1}' AS s
            JOIN (VALUES {2}) AS v(id, m, ketid) ON s.id = v.id)",
            canonical_basis_id, manager->get_path(species, "states"), fmt::join(values, ",")));

        if (result->HasError()) {
            throw cpptrace::runtime_error("Error adding additional kets: " + result->GetError());
        }
    }

    // Ask the table for the extreme values of the quantum numbers
    {
        set_task_status("Validating atomic basis coverage...");

        // Collect the finite quantum-number ranges to validate against the loaded basis. Energy is
        // handled separately because it is compared via its corresponding effective quantum number.
        struct CoverageCheck {
            std::string name;
            std::string column;
            double tolerance{}; // allowed slack between the requested and the available range
            Range<double> range;
        };
        std::vector<CoverageCheck> checks;
        for (const auto &[name, range] : description.quantum_number_ranges) {
            if (!range.is_finite()) {
                continue;
            }
            bool is_expectation_value =
                columns.contains("exp_" + name) && columns.contains("std_" + name);
            std::string column = is_expectation_value ? "exp_" + name : name;
            double tolerance =
                (name == "n" || name == "f" || name == "m" || name == "parity") ? 0.0 : 1.0;
            checks.push_back({name, column, tolerance, range});
        }

        std::string select;
        std::string separator;
        if (description.range_energy.is_finite()) {
            select += separator + "MIN(energy) AS min_energy, MAX(energy) AS max_energy";
            separator = ", ";
        }
        for (const auto &check : checks) {
            select += separator +
                fmt::format("MIN({0}) AS min_{1}, MAX({0}) AS max_{1}", check.column, check.name);
            separator = ", ";
        }

        if (!separator.empty()) {
            auto result =
                con->Query(fmt::format(R"(SELECT {} FROM '{}')", select, canonical_basis_id));

            if (result->HasError()) {
                throw cpptrace::runtime_error("Error querying the database: " + result->GetError());
            }

            auto chunk = result->Fetch();
            const auto &types = result->types;

            for (size_t i = 0; i < chunk->ColumnCount(); i++) {
                if (duckdb::FlatVector::IsNull(chunk->data[i], 0)) {
                    throw std::invalid_argument("No state found.");
                }
            }

            size_t idx = 0;
            if (description.range_energy.is_finite()) {
                auto min_energy = get_entry_as_double(chunk->data[idx], types[idx], 0);
                idx++;
                if (std::sqrt(-1 / (2 * min_energy)) - 1 >
                    std::sqrt(-1 / (2 * description.range_energy.min()))) {
                    SPDLOG_DEBUG("No state found with the requested minimum energy. Requested: {}, "
                                 "found: {}.",
                                 description.range_energy.min(), min_energy);
                }
                auto max_energy = get_entry_as_double(chunk->data[idx], types[idx], 0);
                idx++;
                if (std::sqrt(-1 / (2 * max_energy)) + 1 <
                    std::sqrt(-1 / (2 * description.range_energy.max()))) {
                    SPDLOG_DEBUG("No state found with the requested maximum energy. Requested: {}, "
                                 "found: {}.",
                                 description.range_energy.max(), max_energy);
                }
            }
            for (const auto &check : checks) {
                auto min_value = get_entry_as_double(chunk->data[idx], types[idx], 0);
                idx++;
                if (min_value - check.tolerance > check.range.min()) {
                    SPDLOG_DEBUG("No state found with the requested minimum quantum number {}. "
                                 "Requested: {}, found: {}.",
                                 check.name, check.range.min(), min_value);
                }
                auto max_value = get_entry_as_double(chunk->data[idx], types[idx], 0);
                idx++;
                if (max_value + check.tolerance < check.range.max()) {
                    SPDLOG_DEBUG("No state found with the requested maximum quantum number {}. "
                                 "Requested: {}, found: {}.",
                                 check.name, check.range.max(), max_value);
                }
            }
        }
    }

    // Ask the table for the described states
    set_task_status("Loading atomic basis states...");
    auto result =
        con->Query(fmt::format(R"(SELECT * FROM '{}' ORDER BY ketid ASC)", canonical_basis_id));

    if (result->HasError()) {
        throw cpptrace::runtime_error("Error querying the database: " + result->GetError());
    }

    if (result->RowCount() == 0) {
        throw std::invalid_argument("No state found.");
    }

    // Construct the states. Every column except energy and the raw/encoded id is treated as a
    // quantum number ("id" is the raw states-table id, "ketid" the m-encoded id of the basis
    // state).
    const auto &types = result->types;
    const auto &names = result->names;
    const std::unordered_set<std::string> excluded_columns = {"energy", "id", "ketid"};
    size_t energy_column = get_column_index(names, "energy");
    size_t ketid_column = get_column_index(names, "ketid");

    std::vector<std::shared_ptr<const KetAtom>> kets;
    kets.reserve(result->RowCount());
    double last_energy = std::numeric_limits<double>::lowest();
    double min_quantum_number_nu = std::numeric_limits<double>::max();

    for (auto chunk = result->Fetch(); chunk; chunk = result->Fetch()) {
        set_task_status("Constructing atomic basis...");

        for (size_t i = 0; i < chunk->size(); i++) {
            auto quantum_numbers =
                get_quantum_numbers_from_row(*chunk, types, names, excluded_columns, i);
            double energy =
                get_entry_as_double(chunk->data[energy_column], types[energy_column], i);
            auto id = static_cast<size_t>(
                duckdb::FlatVector::GetData<int64_t>(chunk->data[ketid_column])[i]);

            // Check database consistency
            ensure_consistent_quantum_numbers(quantum_numbers.values.at("f"),
                                              quantum_numbers.values.at("m"));
            if (energy < last_energy) {
                throw std::runtime_error("The states are not sorted by energy.");
            }
            last_energy = energy;

            if (auto it = quantum_numbers.values.find("nu"); it != quantum_numbers.values.end()) {
                min_quantum_number_nu = std::min(min_quantum_number_nu, it->second);
            }

            // Append a new state
            kets.push_back(std::make_shared<const KetAtom>(
                typename KetAtom::Private(), energy, species, std::move(quantum_numbers.values),
                std::move(quantum_numbers.stds), *this, id));
        }
    }

    // Show a warning for low-lying states
    if (min_quantum_number_nu < 25) {
        if (species.ends_with("_mqdt")) {
            SPDLOG_WARN("The multi-channel quantum defect theory might produce inaccurate results "
                        "for effective principal quantum numbers < 25. The models get increasingly "
                        "unreliable for small principal quantum numbers, leading to inaccurate "
                        "matrix elements and energies. Due to missing data, even some states might "
                        "not be present.");
        } else {
            SPDLOG_WARN(
                "The single-channel quantum defect theory can be inaccurate for effective "
                "principal quantum numbers < 25. This can lead to inaccurate matrix elements.");
        }
    }

    return std::make_shared<const BasisAtom<Scalar>>(typename BasisAtom<Scalar>::Private(),
                                                     std::move(kets), std::move(canonical_basis_id),
                                                     *this);
}

template <typename Scalar>
Eigen::SparseMatrix<Scalar, Eigen::RowMajor> Database::get_matrix_elements_in_canonical_basis(
    std::shared_ptr<const BasisAtom<Scalar>> initial_basis,
    std::shared_ptr<const BasisAtom<Scalar>> final_basis, OperatorType type, int q) {
    using real_t = typename traits::NumTraits<Scalar>::real_t;
    using cached_matrix_ptr_t = std::shared_ptr<const cached_matrix_t>;

    if (&initial_basis->get_database() != this || &final_basis->get_database() != this) {
        throw std::invalid_argument(
            "The initial and final bases must belong to the Database instance used for the "
            "matrix element calculation.");
    }
    if (initial_basis->get_species() != final_basis->get_species()) {
        throw std::invalid_argument(fmt::format(
            "The initial and final bases must have the same species, but got '{}' and '{}'.",
            initial_basis->get_species(), final_basis->get_species()));
    }

    std::string specifier;
    int kappa{};
    switch (type) {
    case OperatorType::ELECTRIC_DIPOLE:
        specifier = "matrix_elements_d";
        kappa = 1;
        break;
    case OperatorType::ELECTRIC_QUADRUPOLE:
        specifier = "matrix_elements_q";
        kappa = 2;
        break;
    case OperatorType::ELECTRIC_QUADRUPOLE_ZERO:
        specifier = "matrix_elements_q0";
        kappa = 0;
        break;
    case OperatorType::ELECTRIC_OCTUPOLE:
        specifier = "matrix_elements_o";
        kappa = 3;
        break;
    case OperatorType::MAGNETIC_DIPOLE:
        specifier = "matrix_elements_mu";
        kappa = 1;
        break;
    case OperatorType::ENERGY:
        specifier = "energy";
        kappa = 0;
        break;
    case OperatorType::IDENTITY:
        specifier = "identity";
        kappa = 0;
        break;
    default:
        throw std::invalid_argument("Unknown operator type.");
    }

    std::string canonical_basis_id_initial = initial_basis->get_canonical_basis_id();
    std::string canonical_basis_id_final = final_basis->get_canonical_basis_id();
    std::string cache_key = fmt::format("{}_{}_{}_{}", specifier, q, canonical_basis_id_initial,
                                        canonical_basis_id_final);
    auto &matrix_elements_cache = get_matrix_elements_cache();
    std::promise<cached_matrix_ptr_t> matrix_promise;
    auto [cache_it, inserted] =
        matrix_elements_cache.insert({cache_key, matrix_promise.get_future().share()});

    if (inserted) {
        try {
            Eigen::Index num_rows = final_basis->get_number_of_kets();
            Eigen::Index num_cols = initial_basis->get_number_of_kets();

            std::vector<int> outerIndexPtr;
            std::vector<int> innerIndices;
            std::vector<real_t> values;

            // Check that the specifications are valid
            if (std::abs(q) > kappa) {
                throw std::invalid_argument("Invalid q.");
            }

            // Ask the database for the operator
            set_task_status("Loading matrix elements from database...");
            std::string species = initial_basis->get_species();
            duckdb::unique_ptr<duckdb::MaterializedQueryResult> result;
            if (specifier == "identity") {
                result = con->Query(fmt::format(
                    R"(SELECT s2.ketid AS row, s1.ketid AS col, 1.0::DOUBLE AS val
                    FROM '{}' AS s1
                    INNER JOIN '{}' AS s2 ON s1.ketid = s2.ketid
                    ORDER BY row ASC)",
                    canonical_basis_id_initial, canonical_basis_id_final));
            } else if (specifier == "energy") {
                result = con->Query(fmt::format(
                    R"(SELECT s2.ketid AS row, s1.ketid AS col, s1.energy AS val
                    FROM '{}' AS s1
                    INNER JOIN '{}' AS s2 ON s1.ketid = s2.ketid
                    ORDER BY row ASC)",
                    canonical_basis_id_initial, canonical_basis_id_final));
            } else {
                result = con->Query(fmt::format(
                    R"(WITH s1 AS (
                        SELECT id, f, m, ketid FROM '{}'
                    ),
                    s2 AS (
                        SELECT id, f, m, ketid FROM '{}'
                    ),
                    b AS (
                        SELECT MIN(f) AS min_f, MAX(f) AS max_f,
                        MIN(id) AS min_id, MAX(id) AS max_id
                        FROM (SELECT f, id FROM s1 UNION ALL SELECT f, id FROM s2)
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
                    JOIN s1 ON e.id_initial = s1.id
                    JOIN s2 ON e.id_final = s2.id
                    JOIN w_filtered AS w ON
                    w.f_initial = s1.f AND w.m_initial = s1.m AND
                    w.f_final = s2.f AND w.m_final = s2.m
                    ORDER BY row ASC, col ASC)",
                    canonical_basis_id_initial, canonical_basis_id_final,
                    manager->get_path("misc", "wigner"), kappa, q,
                    manager->get_path(species, specifier)));
            }

            if (result->HasError()) {
                throw cpptrace::runtime_error("Error querying the database: " + result->GetError());
            }

            // Check the types of the columns
            const auto &types = result->types;
            const auto &labels = result->names;
            const std::vector<duckdb::LogicalType> ref_types = {duckdb::LogicalType::BIGINT,
                                                                duckdb::LogicalType::BIGINT,
                                                                duckdb::LogicalType::DOUBLE};
            for (size_t i = 0; i < types.size(); i++) {
                if (types[i] != ref_types[i]) {
                    throw std::runtime_error("Wrong type for '" + labels[i] + "'.");
                }
            }

            set_task_status("Constructing matrix elements...");

            // Construct the matrix
            int num_entries = static_cast<int>(result->RowCount());
            outerIndexPtr.reserve(num_rows + 1);
            innerIndices.reserve(num_entries);
            values.reserve(num_entries);

            int last_row = -1;

            for (auto chunk = result->Fetch(); chunk; chunk = result->Fetch()) {
                auto *chunk_row = duckdb::FlatVector::GetData<int64_t>(chunk->data[0]);
                auto *chunk_col = duckdb::FlatVector::GetData<int64_t>(chunk->data[1]);
                auto *chunk_val = duckdb::FlatVector::GetData<double>(chunk->data[2]);

                for (size_t i = 0; i < chunk->size(); i++) {
                    int row = final_basis->get_ket_index_from_id(chunk_row[i]);
                    if (row != last_row) {
                        if (row < last_row) {
                            throw std::runtime_error("The rows are not sorted.");
                        }
                        for (; last_row < row; last_row++) {
                            outerIndexPtr.push_back(static_cast<int>(innerIndices.size()));
                        }
                    }
                    innerIndices.push_back(initial_basis->get_ket_index_from_id(chunk_col[i]));
                    values.push_back(chunk_val[i]);
                }
            }

            for (; last_row < num_rows + 1; last_row++) {
                outerIndexPtr.push_back(static_cast<int>(innerIndices.size()));
            }

            Eigen::Map<const cached_matrix_t> matrix_map(num_rows, num_cols, values.size(),
                                                         outerIndexPtr.data(), innerIndices.data(),
                                                         values.data());

            auto cached_matrix = std::make_shared<const cached_matrix_t>(matrix_map);
            matrix_promise.set_value(std::move(cached_matrix));
        } catch (...) {
            matrix_promise.set_exception(std::current_exception());
            throw;
        }
    }

    set_task_status("Returning matrix elements in canonical basis...");

    return cache_it->second.get()->template cast<Scalar>();
}

bool Database::get_download_missing() const { return download_missing_; }

bool Database::get_use_cache() const { return use_cache_; }

std::filesystem::path Database::get_database_dir() const { return database_dir_; }

std::string Database::get_versions_info() const { return manager->get_versions_info(); }

Database::matrix_elements_cache_t &Database::get_matrix_elements_cache() {
    static matrix_elements_cache_t matrix_elements_cache;
    return matrix_elements_cache;
}

Database &Database::get_global_instance() {
    return get_global_instance_without_checks(default_download_missing, default_use_cache,
                                              default_database_dir);
}

Database &Database::get_global_instance(bool download_missing) {
    Database &database = get_global_instance_without_checks(download_missing, default_use_cache,
                                                            default_database_dir);
    if (download_missing != database.download_missing_) {
        throw std::invalid_argument(
            "The 'download_missing' argument must not change between calls to the method.");
    }
    return database;
}

Database &Database::get_global_instance(std::filesystem::path database_dir) {
    if (database_dir.empty()) {
        database_dir = default_database_dir;
    }
    Database &database = get_global_instance_without_checks(default_download_missing,
                                                            default_use_cache, database_dir);
    if (!std::filesystem::exists(database_dir) ||
        std::filesystem::canonical(database_dir) != database.database_dir_) {
        throw std::invalid_argument(
            "The 'database_dir' argument must not change between calls to the method.");
    }
    return database;
}

Database &Database::get_global_instance(bool download_missing, bool use_cache,
                                        std::filesystem::path database_dir) {
    if (database_dir.empty()) {
        database_dir = default_database_dir;
    }
    Database &database =
        get_global_instance_without_checks(download_missing, use_cache, database_dir);
    if (download_missing != database.download_missing_ || use_cache != database.use_cache_ ||
        !std::filesystem::exists(database_dir) ||
        std::filesystem::canonical(database_dir) != database.database_dir_) {
        throw std::invalid_argument(
            "The 'download_missing', 'use_cache' and 'database_dir' arguments must not "
            "change between calls to the method.");
    }
    return database;
}

Database &Database::get_global_instance_without_checks(bool download_missing, bool use_cache,
                                                       std::filesystem::path database_dir) {
    static Database database(download_missing, use_cache, std::move(database_dir));
    return database;
}

struct database_dir_noexcept : std::filesystem::path {
    explicit database_dir_noexcept() noexcept try : std
        ::filesystem::path(paths::get_cache_directory() / "database") {}
    catch (...) {
        SPDLOG_ERROR("Error getting the PairInteraction cache directory.");
        std::terminate();
    }
};

const std::filesystem::path Database::default_database_dir = database_dir_noexcept();

// Explicit instantiations
// NOLINTBEGIN(bugprone-macro-parentheses, cppcoreguidelines-macro-usage)
#define INSTANTIATE_GETTERS(SCALAR)                                                                \
    template std::shared_ptr<const BasisAtom<SCALAR>> Database::get_basis<SCALAR>(                 \
        const std::string &species, const AtomDescriptionByRanges &description,                    \
        const std::vector<size_t> &additional_ket_ids);                                            \
    template Eigen::SparseMatrix<SCALAR, Eigen::RowMajor>                                          \
    Database::get_matrix_elements_in_canonical_basis<SCALAR>(                                      \
        std::shared_ptr<const BasisAtom<SCALAR>> initial_basis,                                    \
        std::shared_ptr<const BasisAtom<SCALAR>> final_basis, OperatorType type, int q);
// NOLINTEND(bugprone-macro-parentheses, cppcoreguidelines-macro-usage)

INSTANTIATE_GETTERS(double)
INSTANTIATE_GETTERS(std::complex<double>)

#undef INSTANTIATE_GETTERS
} // namespace pairinteraction
