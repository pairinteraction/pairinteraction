// SPDX-FileCopyrightText: 2024 PairInteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include "pairinteraction/utils/eigen_assertion.hpp"
#include "pairinteraction/utils/traits.hpp"

#include <Eigen/SparseCore>
#include <complex>
#include <cstdint>
#include <filesystem>
#include <future>
#include <memory>
#include <oneapi/tbb.h>
#include <string>
#include <unordered_map>
#include <vector>

namespace duckdb {
class DuckDB;
class Connection;
enum class LogicalTypeId : uint8_t;
} // namespace duckdb

namespace pairinteraction {
enum class OperatorType;

struct AtomDescriptionByParameters;

struct AtomDescriptionByRanges;

class KetAtom;

template <typename Scalar>
class BasisAtom;

class GitHubDownloader;

class ParquetManager;

class Database {
public:
    Database();
    Database(bool download_missing);
    Database(std::filesystem::path database_dir);
    Database(bool download_missing, bool use_cache, std::filesystem::path database_dir);
    ~Database();
    static Database &get_global_instance();
    static Database &get_global_instance(bool download_missing);
    static Database &get_global_instance(std::filesystem::path database_dir);
    static Database &get_global_instance(bool download_missing, bool use_cache,
                                         std::filesystem::path database_dir);

    std::shared_ptr<const KetAtom> get_ket(const std::string &species,
                                           const AtomDescriptionByParameters &description);

    template <typename Scalar>
    std::shared_ptr<const BasisAtom<Scalar>>
    get_basis(const std::string &species, const AtomDescriptionByRanges &description,
              const std::vector<size_t> &additional_ket_ids);

    template <typename Scalar>
    Eigen::SparseMatrix<Scalar, Eigen::RowMajor>
    get_matrix_elements_in_canonical_basis(std::shared_ptr<const BasisAtom<Scalar>> initial_basis,
                                           std::shared_ptr<const BasisAtom<Scalar>> final_basis,
                                           OperatorType type, int q);

    bool get_download_missing() const;
    bool get_use_cache() const;
    std::filesystem::path get_database_dir() const;
    std::string get_versions_info() const;

private:
    struct Table {
        std::filesystem::path local_path{""};
        std::filesystem::path remote_path{""};
        int local_version{-1};
        int remote_version{-1};
    };

    const std::string default_database_repo_host{"https://api.github.com"};
    const std::vector<std::string> default_database_repo_paths{
        "/repos/pairinteraction/database-sqdt/releases/latest",
        "/repos/pairinteraction/database-mqdt/releases/latest"};

    bool download_missing_;
    bool use_cache_;
    std::filesystem::path database_dir_;
    std::unique_ptr<duckdb::DuckDB> db;
    std::unique_ptr<duckdb::Connection> con;
    std::unique_ptr<GitHubDownloader> downloader;
    std::unique_ptr<ParquetManager> manager;

    static constexpr bool default_download_missing{false};
    static constexpr bool default_use_cache{true};
    static const std::filesystem::path default_database_dir;

    using cached_matrix_t = Eigen::SparseMatrix<double, Eigen::RowMajor>;
    using matrix_elements_cache_t = oneapi::tbb::concurrent_unordered_map<
        std::string, std::shared_future<std::shared_ptr<const cached_matrix_t>>>;

    static matrix_elements_cache_t &get_matrix_elements_cache();

    static Database &get_global_instance_without_checks(bool download_missing, bool use_cache,
                                                        std::filesystem::path database_dir);

    void ensure_presence_of_table(const std::string &name);

    // Map from the column names of a states table to the type the column is stored as. The backend
    // never interprets the column names themselves; it only needs to know which columns exist and
    // of which type they are.
    using column_info_t = std::unordered_map<std::string, duckdb::LogicalTypeId>;

    // Return the column info of the states table of the given species, caching the result per table
    // so that the schema is queried from the database at most once.
    const column_info_t &get_states_column_info(const std::string &species);

    oneapi::tbb::concurrent_unordered_map<std::string, column_info_t> column_info_cache;
};

// Extern template declarations
// NOLINTBEGIN(bugprone-macro-parentheses, cppcoreguidelines-macro-usage)
#define EXTERN_GETTERS(SCALAR)                                                                     \
    extern template std::shared_ptr<const BasisAtom<SCALAR>> Database::get_basis<SCALAR>(          \
        const std::string &species, const AtomDescriptionByRanges &description,                    \
        const std::vector<size_t> &additional_ket_ids);                                            \
    extern template Eigen::SparseMatrix<SCALAR, Eigen::RowMajor>                                   \
    Database::get_matrix_elements_in_canonical_basis<SCALAR>(                                      \
        std::shared_ptr<const BasisAtom<SCALAR>> initial_basis,                                    \
        std::shared_ptr<const BasisAtom<SCALAR>> final_basis, OperatorType type, int q);
// NOLINTEND(bugprone-macro-parentheses, cppcoreguidelines-macro-usage)

EXTERN_GETTERS(double)
EXTERN_GETTERS(std::complex<double>)

#undef EXTERN_GETTERS
} // namespace pairinteraction
