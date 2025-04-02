// SPDX-FileCopyrightText: 2024 Pairinteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include "pairinteraction/utils/eigen_assertion.hpp"
#include "pairinteraction/utils/traits.hpp"

#include <Eigen/SparseCore>
#include <complex>
#include <filesystem>
#include <memory>
#include <oneapi/tbb.h>
#include <string>
#include <vector>

namespace duckdb {
class DuckDB;
class Connection;
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
    std::shared_ptr<const BasisAtom<Scalar>> get_basis(const std::string &species,
                                                       const AtomDescriptionByRanges &description,
                                                       std::vector<size_t> additional_ket_ids);

    template <typename Scalar>
    Eigen::SparseMatrix<Scalar, Eigen::RowMajor>
    get_matrix_elements(std::shared_ptr<const BasisAtom<Scalar>> initial_basis,
                        std::shared_ptr<const BasisAtom<Scalar>> final_basis, OperatorType type,
                        int q);

    bool get_download_missing() const;
    bool get_use_cache() const;
    std::filesystem::path get_database_dir() const;

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

    static oneapi::tbb::concurrent_unordered_map<std::string,
                                                 Eigen::SparseMatrix<double, Eigen::RowMajor>> &
    get_matrix_elements_cache();

    static Database &get_global_instance_without_checks(bool download_missing, bool use_cache,
                                                        std::filesystem::path database_dir);

    void ensure_presence_of_table(const std::string &name);
};

// Extern template declarations
// NOLINTBEGIN(bugprone-macro-parentheses, cppcoreguidelines-macro-usage)
#define EXTERN_GETTERS(SCALAR)                                                                     \
    extern template std::shared_ptr<const BasisAtom<SCALAR>> Database::get_basis<SCALAR>(          \
        const std::string &species, const AtomDescriptionByRanges &description,                    \
        std::vector<size_t> additional_ket_ids);                                                   \
    extern template Eigen::SparseMatrix<SCALAR, Eigen::RowMajor>                                   \
    Database::get_matrix_elements<SCALAR>(std::shared_ptr<const BasisAtom<SCALAR>> initial_basis,  \
                                          std::shared_ptr<const BasisAtom<SCALAR>> final_basis,    \
                                          OperatorType type, int q);
// NOLINTEND(bugprone-macro-parentheses, cppcoreguidelines-macro-usage)

EXTERN_GETTERS(double)
EXTERN_GETTERS(std::complex<double>)

#undef EXTERN_GETTERS
} // namespace pairinteraction
