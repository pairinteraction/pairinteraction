#pragma once

#include "pairinteraction/utils/eigen_assertion.hpp"
#include "pairinteraction/utils/traits.hpp"

#include <Eigen/SparseCore>
#include <complex>
#include <filesystem>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

namespace duckdb {
class DuckDB;
class Connection;
} // namespace duckdb

namespace httplib {
class Client;
} // namespace httplib

namespace pairinteraction {
enum class OperatorType;

template <typename Scalar>
struct AtomDescriptionByParameters;

template <typename Scalar>
struct AtomDescriptionByRanges;

template <typename Real>
class KetAtom;

template <typename Scalar>
class BasisAtom;

class Database {
public:
    struct AvailabilitySpecies {
        std::string name;
        bool locally_available;
        bool up_to_date;
        bool fully_downloaded;
    };

    struct AvailabilityWigner {
        bool locally_available;
        bool up_to_date;
    };

    Database();
    Database(bool download_missing);
    Database(std::filesystem::path database_dir);
    Database(bool download_missing, bool wigner_in_memory, std::filesystem::path database_dir);
    ~Database();
    std::vector<AvailabilitySpecies> get_availability_of_species();
    AvailabilityWigner get_availability_of_wigner_table();
    static Database &get_global_instance();
    static Database &get_global_instance(bool download_missing);
    static Database &get_global_instance(std::filesystem::path database_dir);
    static Database &get_global_instance(bool download_missing, bool wigner_in_memory,
                                         std::filesystem::path database_dir);

    template <typename Real>
    std::shared_ptr<const KetAtom<Real>>
    get_ket(std::string species, const AtomDescriptionByParameters<Real> &description);

    template <typename Scalar>
    std::shared_ptr<const BasisAtom<Scalar>> get_basis(
        std::string species,
        const AtomDescriptionByRanges<typename traits::NumTraits<Scalar>::real_t> &description,
        std::vector<size_t> additional_ket_ids);

    template <typename Scalar>
    Eigen::SparseMatrix<Scalar, Eigen::RowMajor>
    get_matrix_elements(std::shared_ptr<const BasisAtom<Scalar>> initial_basis,
                        std::shared_ptr<const BasisAtom<Scalar>> final_basis, OperatorType type,
                        int q);

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

    bool _download_missing;
    bool _wigner_in_memory;
    std::filesystem::path _database_dir;
    std::unique_ptr<duckdb::DuckDB> db;
    std::unique_ptr<duckdb::Connection> con;
    std::vector<httplib::Client> pool;
    std::unordered_map<std::string, Table> tables;

    static constexpr bool default_download_missing{false};
    static constexpr bool default_wigner_in_memory{true};
    static const std::filesystem::path default_database_dir;

    template <typename Real>
    std::unordered_map<std::string, Eigen::SparseMatrix<Real, Eigen::RowMajor>> &
    get_matrix_elements_cache();

    static Database &get_global_instance_without_checks(bool download_missing,
                                                        bool wigner_in_memory,
                                                        std::filesystem::path database_dir);

    void ensure_presence_of_table(const std::string &name);
    void ensure_quantum_number_n_is_allowed(const std::string &name);
};

extern template std::shared_ptr<const KetAtom<float>>
Database::get_ket<float>(std::string species,
                         const AtomDescriptionByParameters<float> &description);
extern template std::shared_ptr<const KetAtom<double>>
Database::get_ket<double>(std::string species,
                          const AtomDescriptionByParameters<double> &description);
extern template std::shared_ptr<const BasisAtom<float>>
Database::get_basis<float>(std::string species, const AtomDescriptionByRanges<float> &description,
                           std::vector<size_t> additional_ket_ids);
extern template std::shared_ptr<const BasisAtom<double>>
Database::get_basis<double>(std::string species, const AtomDescriptionByRanges<double> &description,
                            std::vector<size_t> additional_ket_ids);
extern template std::shared_ptr<const BasisAtom<std::complex<float>>>
Database::get_basis<std::complex<float>>(std::string species,
                                         const AtomDescriptionByRanges<float> &description,
                                         std::vector<size_t> additional_ket_ids);
extern template std::shared_ptr<const BasisAtom<std::complex<double>>>
Database::get_basis<std::complex<double>>(std::string species,
                                          const AtomDescriptionByRanges<double> &description,
                                          std::vector<size_t> additional_ket_ids);
extern template Eigen::SparseMatrix<float, Eigen::RowMajor>
Database::get_matrix_elements<float>(std::shared_ptr<const BasisAtom<float>> initial_basis,
                                     std::shared_ptr<const BasisAtom<float>> final_basis,
                                     OperatorType type, int q);
extern template Eigen::SparseMatrix<double, Eigen::RowMajor>
Database::get_matrix_elements<double>(std::shared_ptr<const BasisAtom<double>> initial_basis,
                                      std::shared_ptr<const BasisAtom<double>> final_basis,
                                      OperatorType type, int q);
extern template Eigen::SparseMatrix<std::complex<float>, Eigen::RowMajor>
Database::get_matrix_elements<std::complex<float>>(
    std::shared_ptr<const BasisAtom<std::complex<float>>> initial_basis,
    std::shared_ptr<const BasisAtom<std::complex<float>>> final_basis, OperatorType type, int q);
extern template Eigen::SparseMatrix<std::complex<double>, Eigen::RowMajor>
Database::get_matrix_elements<std::complex<double>>(
    std::shared_ptr<const BasisAtom<std::complex<double>>> initial_basis,
    std::shared_ptr<const BasisAtom<std::complex<double>>> final_basis, OperatorType type, int q);
} // namespace pairinteraction
