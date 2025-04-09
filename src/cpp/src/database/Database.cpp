// SPDX-FileCopyrightText: 2024 Pairinteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "pairinteraction/database/Database.hpp"

#include "pairinteraction/basis/BasisAtom.hpp"
#include "pairinteraction/database/AtomDescriptionByParameters.hpp"
#include "pairinteraction/database/AtomDescriptionByRanges.hpp"
#include "pairinteraction/database/GitHubDownloader.hpp"
#include "pairinteraction/database/ParquetManager.hpp"
#include "pairinteraction/enums/OperatorType.hpp"
#include "pairinteraction/enums/Parity.hpp"
#include "pairinteraction/ket/KetAtom.hpp"
#include "pairinteraction/utils/hash.hpp"
#include "pairinteraction/utils/id_in_database.hpp"
#include "pairinteraction/utils/paths.hpp"
#include "pairinteraction/utils/streamed.hpp"

#include <cpptrace/cpptrace.hpp>
#include <duckdb.hpp>
#include <fmt/core.h>
#include <fmt/ranges.h>
#include <fstream>
#include <nlohmann/json.hpp>
#include <oneapi/tbb.h>
#include <spdlog/spdlog.h>
#include <system_error>

namespace pairinteraction {
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

std::shared_ptr<const KetAtom> Database::get_ket(const std::string &species,
                                                 const AtomDescriptionByParameters &description) {
    // Check that the specifications are valid
    if (!description.quantum_number_m.has_value()) {
        throw std::invalid_argument("The quantum number m must be specified.");
    }
    if (description.quantum_number_f.has_value() &&
        2 * description.quantum_number_f.value() !=
            std::rint(2 * description.quantum_number_f.value())) {
        throw std::invalid_argument("The quantum number f must be an integer or half-integer.");
    }
    if (description.quantum_number_f.has_value() && description.quantum_number_f.value() < 0) {
        throw std::invalid_argument("The quantum number f must be positive.");
    }
    if (description.quantum_number_j.has_value() &&
        2 * description.quantum_number_j.value() !=
            std::rint(2 * description.quantum_number_j.value())) {
        throw std::invalid_argument("The quantum number j must be an integer or half-integer.");
    }
    if (description.quantum_number_j.has_value() && description.quantum_number_j.value() < 0) {
        throw std::invalid_argument("The quantum number j must be positive.");
    }
    if (description.quantum_number_m.has_value() &&
        2 * description.quantum_number_m.value() !=
            std::rint(2 * description.quantum_number_m.value())) {
        throw std::invalid_argument("The quantum number m must be an integer or half-integer.");
    }

    // Describe the state
    std::string where;
    std::string separator;
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
    if (description.parity != Parity::UNKNOWN) {
        where += separator + fmt::format("parity = {}", fmt::streamed(description.parity));
        separator = " AND ";
    }
    if (description.quantum_number_n.has_value()) {
        where += separator + fmt::format("n = {}", description.quantum_number_n.value());
        separator = " AND ";
    }
    if (description.quantum_number_nu.has_value()) {
        where += separator +
            fmt::format("nu BETWEEN {} AND {}", description.quantum_number_nu.value() - 0.5,
                        description.quantum_number_nu.value() + 0.5);
        separator = " AND ";
    }
    if (description.quantum_number_nui.has_value()) {
        where += separator +
            fmt::format("exp_nui BETWEEN {} AND {}", description.quantum_number_nui.value() - 0.5,
                        description.quantum_number_nui.value() + 0.5);
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
    if (description.quantum_number_l_ryd.has_value()) {
        where += separator +
            fmt::format("exp_l_ryd BETWEEN {} AND {}",
                        description.quantum_number_l_ryd.value() - 0.5,
                        description.quantum_number_l_ryd.value() + 0.5);
        separator = " AND ";
    }
    if (description.quantum_number_j_ryd.has_value()) {
        where += separator +
            fmt::format("exp_j_ryd BETWEEN {} AND {}",
                        description.quantum_number_j_ryd.value() - 0.5,
                        description.quantum_number_j_ryd.value() + 0.5);
        separator = " AND ";
    }
    if (separator.empty()) {
        where += "FALSE";
    }

    std::string orderby;
    separator = "";
    if (description.energy.has_value()) {
        orderby += separator +
            fmt::format("(SQRT(-1/(2*energy)) - {})^2",
                        std::sqrt(-1 / (2 * description.energy.value())));
        separator = " + ";
    }
    if (description.quantum_number_nu.has_value()) {
        orderby += separator + fmt::format("(nu - {})^2", description.quantum_number_nu.value());
        separator = " + ";
    }
    if (description.quantum_number_nui.has_value()) {
        orderby +=
            separator + fmt::format("(exp_nui - {})^2", description.quantum_number_nui.value());
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
    if (description.quantum_number_l_ryd.has_value()) {
        orderby +=
            separator + fmt::format("(exp_l_ryd - {})^2", description.quantum_number_l_ryd.value());
        separator = " + ";
    }
    if (description.quantum_number_j_ryd.has_value()) {
        orderby +=
            separator + fmt::format("(exp_j_ryd - {})^2", description.quantum_number_j_ryd.value());
        separator = " + ";
    }
    if (separator.empty()) {
        orderby += "id";
    }

    // Ask the database for the described state
    auto result = con->Query(fmt::format(
        R"(SELECT energy, f, parity, id, n, nu, exp_nui, std_nui, exp_l, std_l, exp_s, std_s,
        exp_j, std_j, exp_l_ryd, std_l_ryd, exp_j_ryd, std_j_ryd, is_j_total_momentum, is_calculated_with_mqdt, underspecified_channel_contribution, {} AS order_val FROM '{}' WHERE {} ORDER BY order_val ASC LIMIT 2)",
        orderby, manager->get_path(species, "states"), where));

    if (result->HasError()) {
        throw cpptrace::runtime_error("Error querying the database: " + result->GetError());
    }

    if (result->RowCount() == 0) {
        throw std::invalid_argument("No state found.");
    }

    // Check the types of the columns
    const auto &types = result->types;
    const auto &labels = result->names;
    const std::vector<duckdb::LogicalType> ref_types = {
        duckdb::LogicalType::DOUBLE,  duckdb::LogicalType::DOUBLE,  duckdb::LogicalType::BIGINT,
        duckdb::LogicalType::BIGINT,  duckdb::LogicalType::BIGINT,  duckdb::LogicalType::DOUBLE,
        duckdb::LogicalType::DOUBLE,  duckdb::LogicalType::DOUBLE,  duckdb::LogicalType::DOUBLE,
        duckdb::LogicalType::DOUBLE,  duckdb::LogicalType::DOUBLE,  duckdb::LogicalType::DOUBLE,
        duckdb::LogicalType::DOUBLE,  duckdb::LogicalType::DOUBLE,  duckdb::LogicalType::DOUBLE,
        duckdb::LogicalType::DOUBLE,  duckdb::LogicalType::DOUBLE,  duckdb::LogicalType::DOUBLE,
        duckdb::LogicalType::BOOLEAN, duckdb::LogicalType::BOOLEAN, duckdb::LogicalType::DOUBLE,
        duckdb::LogicalType::DOUBLE};

    for (size_t i = 0; i < types.size(); i++) {
        if (types[i] != ref_types[i]) {
            throw std::runtime_error("Wrong type for '" + labels[i] + "'. Got " +
                                     types[i].ToString() + " but expected " +
                                     ref_types[i].ToString());
        }
    }

    // Get the first chunk of the results (the first chunk is sufficient as we need two rows at
    // most)
    auto chunk = result->Fetch();

    // Check that the ket is uniquely specified
    if (chunk->size() > 1) {
        auto order_val_0 = duckdb::FlatVector::GetData<double>(chunk->data[21])[0];
        auto order_val_1 = duckdb::FlatVector::GetData<double>(chunk->data[21])[1];

        if (order_val_1 - order_val_0 <= order_val_0) {
            // Get a list of possible kets
            std::vector<KetAtom> kets;
            kets.reserve(2);
            for (size_t i = 0; i < 2; ++i) {
                auto result_quantum_number_m = description.quantum_number_m.value();
                auto result_energy = duckdb::FlatVector::GetData<double>(chunk->data[0])[i];
                auto result_quantum_number_f =
                    duckdb::FlatVector::GetData<double>(chunk->data[1])[i];
                auto result_parity = duckdb::FlatVector::GetData<int64_t>(chunk->data[2])[i];
                auto result_id = utils::get_linearized_id_in_database(
                    duckdb::FlatVector::GetData<int64_t>(chunk->data[3])[i],
                    result_quantum_number_m);
                auto result_quantum_number_n =
                    duckdb::FlatVector::GetData<int64_t>(chunk->data[4])[i];
                auto result_quantum_number_nu =
                    duckdb::FlatVector::GetData<double>(chunk->data[5])[i];
                auto result_quantum_number_nui_exp =
                    duckdb::FlatVector::GetData<double>(chunk->data[6])[i];
                auto result_quantum_number_nui_std =
                    duckdb::FlatVector::GetData<double>(chunk->data[7])[i];
                auto result_quantum_number_l_exp =
                    duckdb::FlatVector::GetData<double>(chunk->data[8])[i];
                auto result_quantum_number_l_std =
                    duckdb::FlatVector::GetData<double>(chunk->data[9])[i];
                auto result_quantum_number_s_exp =
                    duckdb::FlatVector::GetData<double>(chunk->data[10])[i];
                auto result_quantum_number_s_std =
                    duckdb::FlatVector::GetData<double>(chunk->data[11])[i];
                auto result_quantum_number_j_exp =
                    duckdb::FlatVector::GetData<double>(chunk->data[12])[i];
                auto result_quantum_number_j_std =
                    duckdb::FlatVector::GetData<double>(chunk->data[13])[i];
                auto result_quantum_number_l_ryd_exp =
                    duckdb::FlatVector::GetData<double>(chunk->data[14])[i];
                auto result_quantum_number_l_ryd_std =
                    duckdb::FlatVector::GetData<double>(chunk->data[15])[i];
                auto result_quantum_number_j_ryd_exp =
                    duckdb::FlatVector::GetData<double>(chunk->data[16])[i];
                auto result_quantum_number_j_ryd_std =
                    duckdb::FlatVector::GetData<double>(chunk->data[17])[i];
                auto result_is_j_total_momentum =
                    duckdb::FlatVector::GetData<bool>(chunk->data[18])[i];
                auto result_is_calculated_with_mqdt =
                    duckdb::FlatVector::GetData<bool>(chunk->data[19])[i];
                auto result_underspecified_channel_contribution =
                    duckdb::FlatVector::GetData<double>(chunk->data[20])[i];
                kets.emplace_back(typename KetAtom::Private(), result_energy,
                                  result_quantum_number_f, result_quantum_number_m,
                                  static_cast<Parity>(result_parity), species,
                                  result_quantum_number_n, result_quantum_number_nu,
                                  result_quantum_number_nui_exp, result_quantum_number_nui_std,
                                  result_quantum_number_l_exp, result_quantum_number_l_std,
                                  result_quantum_number_s_exp, result_quantum_number_s_std,
                                  result_quantum_number_j_exp, result_quantum_number_j_std,
                                  result_quantum_number_l_ryd_exp, result_quantum_number_l_ryd_std,
                                  result_quantum_number_j_ryd_exp, result_quantum_number_j_ryd_std,
                                  result_is_j_total_momentum, result_is_calculated_with_mqdt,
                                  result_underspecified_channel_contribution, *this, result_id);
            }

            // Throw an error with the possible kets
            throw std::invalid_argument(
                fmt::format("The ket is not uniquely specified. Possible kets are:\n{}\n{}",
                            fmt::streamed(kets[0]), fmt::streamed(kets[1])));
        }
    }

    // Construct the state
    auto result_quantum_number_m = description.quantum_number_m.value();
    auto result_energy = duckdb::FlatVector::GetData<double>(chunk->data[0])[0];
    auto result_quantum_number_f = duckdb::FlatVector::GetData<double>(chunk->data[1])[0];
    auto result_parity = duckdb::FlatVector::GetData<int64_t>(chunk->data[2])[0];
    auto result_id = utils::get_linearized_id_in_database(
        duckdb::FlatVector::GetData<int64_t>(chunk->data[3])[0], result_quantum_number_m);
    auto result_quantum_number_n = duckdb::FlatVector::GetData<int64_t>(chunk->data[4])[0];
    auto result_quantum_number_nu = duckdb::FlatVector::GetData<double>(chunk->data[5])[0];
    auto result_quantum_number_nui_exp = duckdb::FlatVector::GetData<double>(chunk->data[6])[0];
    auto result_quantum_number_nui_std = duckdb::FlatVector::GetData<double>(chunk->data[7])[0];
    auto result_quantum_number_l_exp = duckdb::FlatVector::GetData<double>(chunk->data[8])[0];
    auto result_quantum_number_l_std = duckdb::FlatVector::GetData<double>(chunk->data[9])[0];
    auto result_quantum_number_s_exp = duckdb::FlatVector::GetData<double>(chunk->data[10])[0];
    auto result_quantum_number_s_std = duckdb::FlatVector::GetData<double>(chunk->data[11])[0];
    auto result_quantum_number_j_exp = duckdb::FlatVector::GetData<double>(chunk->data[12])[0];
    auto result_quantum_number_j_std = duckdb::FlatVector::GetData<double>(chunk->data[13])[0];
    auto result_quantum_number_l_ryd_exp = duckdb::FlatVector::GetData<double>(chunk->data[14])[0];
    auto result_quantum_number_l_ryd_std = duckdb::FlatVector::GetData<double>(chunk->data[15])[0];
    auto result_quantum_number_j_ryd_exp = duckdb::FlatVector::GetData<double>(chunk->data[16])[0];
    auto result_quantum_number_j_ryd_std = duckdb::FlatVector::GetData<double>(chunk->data[17])[0];
    auto result_is_j_total_momentum = duckdb::FlatVector::GetData<bool>(chunk->data[18])[0];
    auto result_is_calculated_with_mqdt = duckdb::FlatVector::GetData<bool>(chunk->data[19])[0];
    auto result_underspecified_channel_contribution =
        duckdb::FlatVector::GetData<double>(chunk->data[20])[0];

    // Check the quantum number m
    if (std::abs(result_quantum_number_m) > result_quantum_number_f) {
        throw std::invalid_argument(
            "The absolute value of the quantum number m must be less than or equal to f.");
    }
    if (result_quantum_number_f + result_quantum_number_m !=
        std::rint(result_quantum_number_f + result_quantum_number_m)) {
        throw std::invalid_argument(
            "The quantum numbers f and m must be both either integers or half-integers.");
    }

#ifndef NDEBUG
    // Check database consistency
    if (result_is_j_total_momentum && result_quantum_number_f != result_quantum_number_j_exp) {
        throw std::runtime_error("If j is the total momentum, f must be equal to j.");
    }
#endif

    return std::make_shared<const KetAtom>(
        typename KetAtom::Private(), result_energy, result_quantum_number_f,
        result_quantum_number_m, static_cast<Parity>(result_parity), species,
        result_quantum_number_n, result_quantum_number_nu, result_quantum_number_nui_exp,
        result_quantum_number_nui_std, result_quantum_number_l_exp, result_quantum_number_l_std,
        result_quantum_number_s_exp, result_quantum_number_s_std, result_quantum_number_j_exp,
        result_quantum_number_j_std, result_quantum_number_l_ryd_exp,
        result_quantum_number_l_ryd_std, result_quantum_number_j_ryd_exp,
        result_quantum_number_j_ryd_std, result_is_j_total_momentum, result_is_calculated_with_mqdt,
        result_underspecified_channel_contribution, *this, result_id);
}

template <typename Scalar>
std::shared_ptr<const BasisAtom<Scalar>>
Database::get_basis(const std::string &species, const AtomDescriptionByRanges &description,
                    std::vector<size_t> additional_ket_ids) {
    // Describe the states
    std::string where = "(";
    std::string separator;
    if (description.parity != Parity::UNKNOWN) {
        where += separator + fmt::format("parity = {}", fmt::streamed(description.parity));
        separator = " AND ";
    }
    if (description.range_energy.is_finite()) {
        where += separator +
            fmt::format("energy BETWEEN {} AND {}", description.range_energy.min(),
                        description.range_energy.max());
        separator = " AND ";
    }
    if (description.range_quantum_number_f.is_finite()) {
        where += separator +
            fmt::format("f BETWEEN {} AND {}", description.range_quantum_number_f.min(),
                        description.range_quantum_number_f.max());
        separator = " AND ";
    }
    if (description.range_quantum_number_m.is_finite()) {
        where += separator +
            fmt::format("m BETWEEN {} AND {}", description.range_quantum_number_m.min(),
                        description.range_quantum_number_m.max());
        separator = " AND ";
    }
    if (description.range_quantum_number_n.is_finite()) {
        where += separator +
            fmt::format("n BETWEEN {} AND {}", description.range_quantum_number_n.min(),
                        description.range_quantum_number_n.max());
        separator = " AND ";
    }
    if (description.range_quantum_number_nu.is_finite()) {
        where += separator +
            fmt::format("nu BETWEEN {} AND {}", description.range_quantum_number_nu.min(),
                        description.range_quantum_number_nu.max());
        separator = " AND ";
    }
    if (description.range_quantum_number_nui.is_finite()) {
        where += separator +
            fmt::format("exp_nui BETWEEN {}-2*std_nui AND {}+2*std_nui",
                        description.range_quantum_number_nui.min(),
                        description.range_quantum_number_nui.max());
        separator = " AND ";
    }
    if (description.range_quantum_number_l.is_finite()) {
        where += separator +
            fmt::format("exp_l BETWEEN {}-2*std_l AND {}+2*std_l",
                        description.range_quantum_number_l.min(),
                        description.range_quantum_number_l.max());
        separator = " AND ";
    }
    if (description.range_quantum_number_s.is_finite()) {
        where += separator +
            fmt::format("exp_s BETWEEN {}-2*std_s AND {}+2*std_s",
                        description.range_quantum_number_s.min(),
                        description.range_quantum_number_s.max());
        separator = " AND ";
    }
    if (description.range_quantum_number_j.is_finite()) {
        where += separator +
            fmt::format("exp_j BETWEEN {}-2*std_j AND {}+2*std_j",
                        description.range_quantum_number_j.min(),
                        description.range_quantum_number_j.max());
        separator = " AND ";
    }
    if (description.range_quantum_number_l_ryd.is_finite()) {
        where += separator +
            fmt::format("exp_l_ryd BETWEEN {}-2*std_l_ryd AND {}+2*std_l_ryd",
                        description.range_quantum_number_l_ryd.min(),
                        description.range_quantum_number_l_ryd.max());
        separator = " AND ";
    }
    if (description.range_quantum_number_j_ryd.is_finite()) {
        where += separator +
            fmt::format("exp_j_ryd BETWEEN {}-2*std_j_ryd AND {}+2*std_j_ryd",
                        description.range_quantum_number_j_ryd.min(),
                        description.range_quantum_number_j_ryd.max());
        separator = " AND ";
    }
    if (separator.empty()) {
        where += "FALSE";
    }
    where += ")";
    if (!additional_ket_ids.empty()) {
        where += fmt::format(" OR {} IN ({})", utils::SQL_TERM_FOR_LINEARIZED_ID_IN_DATABASE,
                             fmt::join(additional_ket_ids, ","));
    }

    // Create a table containing the described states
    std::string id_of_kets;
    {
        auto result = con->Query(R"(SELECT UUID()::varchar)");
        if (result->HasError()) {
            throw cpptrace::runtime_error("Error selecting id_of_kets: " + result->GetError());
        }
        id_of_kets =
            duckdb::FlatVector::GetData<duckdb::string_t>(result->Fetch()->data[0])[0].GetString();
    }
    {
        auto result = con->Query(fmt::format(
            R"(CREATE TEMP TABLE '{}' AS SELECT *, {} AS ketid FROM (
                SELECT *,
                UNNEST(list_transform(generate_series(0,(2*f)::bigint),
                x -> x::double-f)) AS m FROM '{}'
            ) WHERE {})",
            id_of_kets, utils::SQL_TERM_FOR_LINEARIZED_ID_IN_DATABASE,
            manager->get_path(species, "states"), where));

        if (result->HasError()) {
            throw cpptrace::runtime_error("Error creating table: " + result->GetError());
        }
    }

    // Ask the table for the extreme values of the quantum numbers
    {
        std::string select;
        std::string separator;
        if (description.range_energy.is_finite()) {
            select += separator + "MIN(energy) AS min_energy, MAX(energy) AS max_energy";
            separator = ", ";
        }
        if (description.range_quantum_number_f.is_finite()) {
            select += separator + "MIN(f) AS min_f, MAX(f) AS max_f";
            separator = ", ";
        }
        if (description.range_quantum_number_m.is_finite()) {
            select += separator + "MIN(m) AS min_m, MAX(m) AS max_m";
            separator = ", ";
        }
        if (description.range_quantum_number_n.is_finite()) {
            select += separator + "MIN(n) AS min_n, MAX(n) AS max_n";
            separator = ", ";
        }
        if (description.range_quantum_number_nu.is_finite()) {
            select += separator + "MIN(nu) AS min_nu, MAX(nu) AS max_nu";
            separator = ", ";
        }
        if (description.range_quantum_number_nui.is_finite()) {
            select += separator + "MIN(exp_nui) AS min_nui, MAX(exp_nui) AS max_nui";
            separator = ", ";
        }
        if (description.range_quantum_number_l.is_finite()) {
            select += separator + "MIN(exp_l) AS min_l, MAX(exp_l) AS max_l";
            separator = ", ";
        }
        if (description.range_quantum_number_s.is_finite()) {
            select += separator + "MIN(exp_s) AS min_s, MAX(exp_s) AS max_s";
            separator = ", ";
        }
        if (description.range_quantum_number_j.is_finite()) {
            select += separator + "MIN(exp_j) AS min_j, MAX(exp_j) AS max_j";
            separator = ", ";
        }
        if (description.range_quantum_number_l_ryd.is_finite()) {
            select += separator + "MIN(exp_l_ryd) AS min_l_ryd, MAX(exp_l_ryd) AS max_l_ryd";
            separator = ", ";
        }
        if (description.range_quantum_number_j_ryd.is_finite()) {
            select += separator + "MIN(exp_j_ryd) AS min_j_ryd, MAX(exp_j_ryd) AS max_j_ryd";
            separator = ", ";
        }

        if (!separator.empty()) {
            auto result = con->Query(fmt::format(R"(SELECT {} FROM '{}')", select, id_of_kets));

            if (result->HasError()) {
                throw cpptrace::runtime_error("Error querying the database: " + result->GetError());
            }

            auto chunk = result->Fetch();

            for (size_t i = 0; i < chunk->ColumnCount(); i++) {
                if (duckdb::FlatVector::IsNull(chunk->data[i], 0)) {
                    throw std::invalid_argument("No state found.");
                }
            }

            size_t idx = 0;
            if (description.range_energy.is_finite()) {
                auto min_energy = duckdb::FlatVector::GetData<double>(chunk->data[idx++])[0];
                if (std::sqrt(-1 / (2 * min_energy)) - 1 >
                    std::sqrt(-1 / (2 * description.range_energy.min()))) {
                    SPDLOG_DEBUG("No state found with the requested minimum energy. Requested: {}, "
                                 "found: {}.",
                                 description.range_energy.min(), min_energy);
                }
                auto max_energy = duckdb::FlatVector::GetData<double>(chunk->data[idx++])[0];
                if (std::sqrt(-1 / (2 * max_energy)) + 1 <
                    std::sqrt(-1 / (2 * description.range_energy.max()))) {
                    SPDLOG_DEBUG("No state found with the requested maximum energy. Requested: {}, "
                                 "found: {}.",
                                 description.range_energy.max(), max_energy);
                }
            }
            if (description.range_quantum_number_f.is_finite()) {
                auto min_f = duckdb::FlatVector::GetData<double>(chunk->data[idx++])[0];
                if (min_f > description.range_quantum_number_f.min()) {
                    SPDLOG_DEBUG("No state found with the requested minimum quantum number f. "
                                 "Requested: {}, found: {}.",
                                 description.range_quantum_number_f.min(), min_f);
                }
                auto max_f = duckdb::FlatVector::GetData<double>(chunk->data[idx++])[0];
                if (max_f < description.range_quantum_number_f.max()) {
                    SPDLOG_DEBUG("No state found with the requested maximum quantum number f. "
                                 "Requested: {}, found: {}.",
                                 description.range_quantum_number_f.max(), max_f);
                }
            }
            if (description.range_quantum_number_m.is_finite()) {
                auto min_m = duckdb::FlatVector::GetData<double>(chunk->data[idx++])[0];
                if (min_m > description.range_quantum_number_m.min()) {
                    SPDLOG_DEBUG("No state found with the requested minimum quantum number m. "
                                 "Requested: {}, found: {}.",
                                 description.range_quantum_number_m.min(), min_m);
                }
                auto max_m = duckdb::FlatVector::GetData<double>(chunk->data[idx++])[0];
                if (max_m < description.range_quantum_number_m.max()) {
                    SPDLOG_DEBUG("No state found with the requested maximum quantum number m. "
                                 "Requested: {}, found: {}.",
                                 description.range_quantum_number_m.max(), max_m);
                }
            }
            if (description.range_quantum_number_n.is_finite()) {
                auto min_n = duckdb::FlatVector::GetData<int64_t>(chunk->data[idx++])[0];
                if (min_n > description.range_quantum_number_n.min()) {
                    SPDLOG_DEBUG("No state found with the requested minimum quantum number n. "
                                 "Requested: {}, found: {}.",
                                 description.range_quantum_number_n.min(), min_n);
                }
                auto max_n = duckdb::FlatVector::GetData<int64_t>(chunk->data[idx++])[0];
                if (max_n < description.range_quantum_number_n.max()) {
                    SPDLOG_DEBUG("No state found with the requested maximum quantum number n. "
                                 "Requested: {}, found: {}.",
                                 description.range_quantum_number_n.max(), max_n);
                }
            }
            if (description.range_quantum_number_nu.is_finite()) {
                auto min_nu = duckdb::FlatVector::GetData<double>(chunk->data[idx++])[0];
                if (min_nu - 1 > description.range_quantum_number_nu.min()) {
                    SPDLOG_DEBUG("No state found with the requested minimum quantum number nu. "
                                 "Requested: {}, found: {}.",
                                 description.range_quantum_number_nu.min(), min_nu);
                }
                auto max_nu = duckdb::FlatVector::GetData<double>(chunk->data[idx++])[0];
                if (max_nu + 1 < description.range_quantum_number_nu.max()) {
                    SPDLOG_DEBUG("No state found with the requested maximum quantum number nu. "
                                 "Requested: {}, found: {}.",
                                 description.range_quantum_number_nu.max(), max_nu);
                }
            }
            if (description.range_quantum_number_nui.is_finite()) {
                auto min_nui = duckdb::FlatVector::GetData<double>(chunk->data[idx++])[0];
                if (min_nui - 1 > description.range_quantum_number_nui.min()) {
                    SPDLOG_DEBUG("No state found with the requested minimum quantum number nui. "
                                 "Requested: {}, found: {}.",
                                 description.range_quantum_number_nui.min(), min_nui);
                }
                auto max_nui = duckdb::FlatVector::GetData<double>(chunk->data[idx++])[0];
                if (max_nui + 1 < description.range_quantum_number_nui.max()) {
                    SPDLOG_DEBUG("No state found with the requested maximum quantum number nui. "
                                 "Requested: {}, found: {}.",
                                 description.range_quantum_number_nui.max(), max_nui);
                }
            }
            if (description.range_quantum_number_l.is_finite()) {
                auto min_l = duckdb::FlatVector::GetData<double>(chunk->data[idx++])[0];
                if (min_l - 1 > description.range_quantum_number_l.min()) {
                    SPDLOG_DEBUG("No state found with the requested minimum quantum number l. "
                                 "Requested: {}, found: {}.",
                                 description.range_quantum_number_l.min(), min_l);
                }
                auto max_l = duckdb::FlatVector::GetData<double>(chunk->data[idx++])[0];
                if (max_l + 1 < description.range_quantum_number_l.max()) {
                    SPDLOG_DEBUG("No state found with the requested maximum quantum number l. "
                                 "Requested: {}, found: {}.",
                                 description.range_quantum_number_l.max(), max_l);
                }
            }
            if (description.range_quantum_number_s.is_finite()) {
                auto min_s = duckdb::FlatVector::GetData<double>(chunk->data[idx++])[0];
                if (min_s - 1 > description.range_quantum_number_s.min()) {
                    SPDLOG_DEBUG("No state found with the requested minimum quantum number s. "
                                 "Requested: {}, found: {}.",
                                 description.range_quantum_number_s.min(), min_s);
                }
                auto max_s = duckdb::FlatVector::GetData<double>(chunk->data[idx++])[0];
                if (max_s + 1 < description.range_quantum_number_s.max()) {
                    SPDLOG_DEBUG("No state found with the requested maximum quantum number s. "
                                 "Requested: {}, found: {}.",
                                 description.range_quantum_number_s.max(), max_s);
                }
            }
            if (description.range_quantum_number_j.is_finite()) {
                auto min_j = duckdb::FlatVector::GetData<double>(chunk->data[idx++])[0];
                if (min_j - 1 > description.range_quantum_number_j.min()) {
                    SPDLOG_DEBUG("No state found with the requested minimum quantum number j. "
                                 "Requested: {}, found: {}.",
                                 description.range_quantum_number_j.min(), min_j);
                }
                auto max_j = duckdb::FlatVector::GetData<double>(chunk->data[idx++])[0];
                if (max_j + 1 < description.range_quantum_number_j.max()) {
                    SPDLOG_DEBUG("No state found with the requested maximum quantum number j. "
                                 "Requested: {}, found: {}.",
                                 description.range_quantum_number_j.max(), max_j);
                }
            }
            if (description.range_quantum_number_l_ryd.is_finite()) {
                auto min_l_ryd = duckdb::FlatVector::GetData<double>(chunk->data[idx++])[0];
                if (min_l_ryd - 1 > description.range_quantum_number_l_ryd.min()) {
                    SPDLOG_DEBUG("No state found with the requested minimum quantum number l_ryd. "
                                 "Requested: {}, found: {}.",
                                 description.range_quantum_number_l_ryd.min(), min_l_ryd);
                }
                auto max_l_ryd = duckdb::FlatVector::GetData<double>(chunk->data[idx++])[0];
                if (max_l_ryd + 1 < description.range_quantum_number_l_ryd.max()) {
                    SPDLOG_DEBUG("No state found with the requested maximum quantum number l_ryd. "
                                 "Requested: {}, found: {}.",
                                 description.range_quantum_number_l_ryd.max(), max_l_ryd);
                }
            }
            if (description.range_quantum_number_j_ryd.is_finite()) {
                auto min_j_ryd = duckdb::FlatVector::GetData<double>(chunk->data[idx++])[0];
                if (min_j_ryd - 1 > description.range_quantum_number_j_ryd.min()) {
                    SPDLOG_DEBUG("No state found with the requested minimum quantum number j_ryd. "
                                 "Requested: {}, found: {}.",
                                 description.range_quantum_number_j_ryd.min(), min_j_ryd);
                }
                auto max_j_ryd = duckdb::FlatVector::GetData<double>(chunk->data[idx++])[0];
                if (max_j_ryd + 1 < description.range_quantum_number_j_ryd.max()) {
                    SPDLOG_DEBUG("No state found with the requested maximum quantum number j_ryd. "
                                 "Requested: {}, found: {}.",
                                 description.range_quantum_number_j_ryd.max(), max_j_ryd);
                }
            }
        }
    }

    // Ask the table for the described states
    auto result = con->Query(fmt::format(
        R"(SELECT energy, f, m, parity, ketid, n, nu, exp_nui, std_nui, exp_l, std_l,
        exp_s, std_s, exp_j, std_j, exp_l_ryd, std_l_ryd, exp_j_ryd, std_j_ryd, is_j_total_momentum, is_calculated_with_mqdt, underspecified_channel_contribution FROM '{}' ORDER BY ketid ASC)",
        id_of_kets));

    if (result->HasError()) {
        throw cpptrace::runtime_error("Error querying the database: " + result->GetError());
    }

    if (result->RowCount() == 0) {
        throw std::invalid_argument("No state found.");
    }

    // Check the types of the columns
    const auto &types = result->types;
    const auto &labels = result->names;
    const std::vector<duckdb::LogicalType> ref_types = {
        duckdb::LogicalType::DOUBLE, duckdb::LogicalType::DOUBLE,  duckdb::LogicalType::DOUBLE,
        duckdb::LogicalType::BIGINT, duckdb::LogicalType::BIGINT,  duckdb::LogicalType::BIGINT,
        duckdb::LogicalType::DOUBLE, duckdb::LogicalType::DOUBLE,  duckdb::LogicalType::DOUBLE,
        duckdb::LogicalType::DOUBLE, duckdb::LogicalType::DOUBLE,  duckdb::LogicalType::DOUBLE,
        duckdb::LogicalType::DOUBLE, duckdb::LogicalType::DOUBLE,  duckdb::LogicalType::DOUBLE,
        duckdb::LogicalType::DOUBLE, duckdb::LogicalType::DOUBLE,  duckdb::LogicalType::DOUBLE,
        duckdb::LogicalType::DOUBLE, duckdb::LogicalType::BOOLEAN, duckdb::LogicalType::BOOLEAN,
        duckdb::LogicalType::DOUBLE};

    for (size_t i = 0; i < types.size(); i++) {
        if (types[i] != ref_types[i]) {
            throw std::runtime_error("Wrong type for '" + labels[i] + "'. Got " +
                                     types[i].ToString() + " but expected " +
                                     ref_types[i].ToString());
        }
    }

    // Construct the states
    std::vector<std::shared_ptr<const KetAtom>> kets;
    kets.reserve(result->RowCount());
#ifndef NDEBUG
    double last_energy = std::numeric_limits<double>::lowest();
#endif

    for (auto chunk = result->Fetch(); chunk; chunk = result->Fetch()) {

        auto *chunk_energy = duckdb::FlatVector::GetData<double>(chunk->data[0]);
        auto *chunk_quantum_number_f = duckdb::FlatVector::GetData<double>(chunk->data[1]);
        auto *chunk_quantum_number_m = duckdb::FlatVector::GetData<double>(chunk->data[2]);
        auto *chunk_parity = duckdb::FlatVector::GetData<int64_t>(chunk->data[3]);
        auto *chunk_id = duckdb::FlatVector::GetData<int64_t>(chunk->data[4]);
        auto *chunk_quantum_number_n = duckdb::FlatVector::GetData<int64_t>(chunk->data[5]);
        auto *chunk_quantum_number_nu = duckdb::FlatVector::GetData<double>(chunk->data[6]);
        auto *chunk_quantum_number_nui_exp = duckdb::FlatVector::GetData<double>(chunk->data[7]);
        auto *chunk_quantum_number_nui_std = duckdb::FlatVector::GetData<double>(chunk->data[8]);
        auto *chunk_quantum_number_l_exp = duckdb::FlatVector::GetData<double>(chunk->data[9]);
        auto *chunk_quantum_number_l_std = duckdb::FlatVector::GetData<double>(chunk->data[10]);
        auto *chunk_quantum_number_s_exp = duckdb::FlatVector::GetData<double>(chunk->data[11]);
        auto *chunk_quantum_number_s_std = duckdb::FlatVector::GetData<double>(chunk->data[12]);
        auto *chunk_quantum_number_j_exp = duckdb::FlatVector::GetData<double>(chunk->data[13]);
        auto *chunk_quantum_number_j_std = duckdb::FlatVector::GetData<double>(chunk->data[14]);
        auto *chunk_quantum_number_l_ryd_exp = duckdb::FlatVector::GetData<double>(chunk->data[15]);
        auto *chunk_quantum_number_l_ryd_std = duckdb::FlatVector::GetData<double>(chunk->data[16]);
        auto *chunk_quantum_number_j_ryd_exp = duckdb::FlatVector::GetData<double>(chunk->data[17]);
        auto *chunk_quantum_number_j_ryd_std = duckdb::FlatVector::GetData<double>(chunk->data[18]);
        auto *chunk_is_j_total_momentum = duckdb::FlatVector::GetData<bool>(chunk->data[19]);
        auto *chunk_is_calculated_with_mqdt = duckdb::FlatVector::GetData<bool>(chunk->data[20]);
        auto *chunk_underspecified_channel_contribution =
            duckdb::FlatVector::GetData<double>(chunk->data[21]);

        for (size_t i = 0; i < chunk->size(); i++) {

#ifndef NDEBUG
            // Check database consistency
            if (chunk_is_j_total_momentum[i] &&
                chunk_quantum_number_f[i] != chunk_quantum_number_j_exp[i]) {
                throw std::runtime_error("If j is the total momentum, f must be equal to j.");
            }
            if (chunk_energy[i] < last_energy) {
                throw std::runtime_error("The states are not sorted by energy.");
            }
            last_energy = chunk_energy[i];
#endif

            // Append a new state
            kets.push_back(std::make_shared<const KetAtom>(
                typename KetAtom::Private(), chunk_energy[i], chunk_quantum_number_f[i],
                chunk_quantum_number_m[i], static_cast<Parity>(chunk_parity[i]), species,
                chunk_quantum_number_n[i], chunk_quantum_number_nu[i],
                chunk_quantum_number_nui_exp[i], chunk_quantum_number_nui_std[i],
                chunk_quantum_number_l_exp[i], chunk_quantum_number_l_std[i],
                chunk_quantum_number_s_exp[i], chunk_quantum_number_s_std[i],
                chunk_quantum_number_j_exp[i], chunk_quantum_number_j_std[i],
                chunk_quantum_number_l_ryd_exp[i], chunk_quantum_number_l_ryd_std[i],
                chunk_quantum_number_j_ryd_exp[i], chunk_quantum_number_j_ryd_std[i],
                chunk_is_j_total_momentum[i], chunk_is_calculated_with_mqdt[i],
                chunk_underspecified_channel_contribution[i], *this, chunk_id[i]));
        }
    }

    return std::make_shared<const BasisAtom<Scalar>>(typename BasisAtom<Scalar>::Private(),
                                                     std::move(kets), std::move(id_of_kets), *this);
}

template <typename Scalar>
Eigen::SparseMatrix<Scalar, Eigen::RowMajor>
Database::get_matrix_elements(std::shared_ptr<const BasisAtom<Scalar>> initial_basis,
                              std::shared_ptr<const BasisAtom<Scalar>> final_basis,
                              OperatorType type, int q) {
    using real_t = typename traits::NumTraits<Scalar>::real_t;

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

    if (initial_basis->get_id_of_kets() != final_basis->get_id_of_kets()) {
        throw std::invalid_argument(
            "The initial and final basis must be expressed using the same kets.");
    }
    std::string id_of_kets = initial_basis->get_id_of_kets();
    std::string cache_key = fmt::format("{}_{}_{}", specifier, q, id_of_kets);

    if (!get_matrix_elements_cache().contains(cache_key)) {
        Eigen::Index dim = initial_basis->get_number_of_kets();

        std::vector<int> outerIndexPtr;
        std::vector<int> innerIndices;
        std::vector<real_t> values;

        if (specifier == "identity") {
            outerIndexPtr.reserve(dim + 1);
            innerIndices.reserve(dim);
            values.reserve(dim);

            for (int i = 0; i < dim; i++) {
                outerIndexPtr.push_back(static_cast<int>(innerIndices.size()));
                innerIndices.push_back(i);
                values.push_back(1);
            }
            outerIndexPtr.push_back(static_cast<int>(innerIndices.size()));

        } else {
            // Check that the specifications are valid
            if (std::abs(q) > kappa) {
                throw std::invalid_argument("Invalid q.");
            }

            // Ask the database for the operator
            std::string species = initial_basis->get_species();
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
                    id_of_kets, manager->get_path("misc", "wigner"), kappa, q,
                    manager->get_path(species, specifier)));
            } else {
                result = con->Query(fmt::format(
                    R"(SELECT ketid as row, ketid as col, energy as val FROM '{}' ORDER BY row ASC)",
                    id_of_kets));
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

            // Construct the matrix
            int num_entries = static_cast<int>(result->RowCount());
            outerIndexPtr.reserve(dim + 1);
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

            for (; last_row < dim + 1; last_row++) {
                outerIndexPtr.push_back(static_cast<int>(innerIndices.size()));
            }
        }

        Eigen::Map<const Eigen::SparseMatrix<real_t, Eigen::RowMajor>> matrix_map(
            dim, dim, values.size(), outerIndexPtr.data(), innerIndices.data(), values.data());

        // Cache the matrix
        get_matrix_elements_cache()[cache_key] = matrix_map;
    }

    // Construct the operator and return it
    return final_basis->get_coefficients().adjoint() *
        get_matrix_elements_cache()[cache_key].template cast<Scalar>() *
        initial_basis->get_coefficients();
}

bool Database::get_download_missing() const { return download_missing_; }

bool Database::get_use_cache() const { return use_cache_; }

std::filesystem::path Database::get_database_dir() const { return database_dir_; }

oneapi::tbb::concurrent_unordered_map<std::string, Eigen::SparseMatrix<double, Eigen::RowMajor>> &
Database::get_matrix_elements_cache() {
    static oneapi::tbb::concurrent_unordered_map<std::string,
                                                 Eigen::SparseMatrix<double, Eigen::RowMajor>>
        matrix_elements_cache;
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
        SPDLOG_ERROR("Error getting the pairinteraction cache directory.");
        std::terminate();
    }
};

const std::filesystem::path Database::default_database_dir = database_dir_noexcept();

// Explicit instantiations
// NOLINTBEGIN(bugprone-macro-parentheses, cppcoreguidelines-macro-usage)
#define INSTANTIATE_GETTERS(SCALAR)                                                                \
    template std::shared_ptr<const BasisAtom<SCALAR>> Database::get_basis<SCALAR>(                 \
        const std::string &species, const AtomDescriptionByRanges &description,                    \
        std::vector<size_t> additional_ket_ids);                                                   \
    template Eigen::SparseMatrix<SCALAR, Eigen::RowMajor> Database::get_matrix_elements<SCALAR>(   \
        std::shared_ptr<const BasisAtom<SCALAR>> initial_basis,                                    \
        std::shared_ptr<const BasisAtom<SCALAR>> final_basis, OperatorType type, int q);
// NOLINTEND(bugprone-macro-parentheses, cppcoreguidelines-macro-usage)

INSTANTIATE_GETTERS(double)
INSTANTIATE_GETTERS(std::complex<double>)

#undef INSTANTIATE_GETTERS
} // namespace pairinteraction
