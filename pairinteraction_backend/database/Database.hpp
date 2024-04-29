#pragma once

#include <Eigen/SparseCore>
#include <complex>
#include <filesystem>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

#include "enums/OperatorType.hpp"

template <typename Scalar>
struct AtomDescriptionByParameters;

template <typename Scalar>
struct AtomDescriptionByRanges;

template <typename Real>
class KetAtom;

template <typename Real>
class BasisAtom;

namespace duckdb {
class DuckDB;
class Connection;
} // namespace duckdb

namespace httplib {
class Client;
} // namespace httplib

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

    Database(bool auto_update = true);
    ~Database();
    std::vector<AvailabilitySpecies> get_availability_of_species();
    AvailabilityWigner get_availability_of_wigner_table();
    static Database &get_global_instance();

    template <typename Real>
    KetAtom<Real> get_ket(std::string species,
                          const AtomDescriptionByParameters<Real> &description);

    template <typename Scalar>
    BasisAtom<Scalar> get_basis(std::string species,
                                const AtomDescriptionByRanges<Scalar> &description,
                                std::vector<size_t> additional_ket_ids);

    template <typename Scalar>
    Eigen::SparseMatrix<Scalar, Eigen::RowMajor> get_operator(const BasisAtom<Scalar> &basis,
                                                              OperatorType type, int q);

private:
    struct Table {
        std::filesystem::path local_path{""};
        std::filesystem::path remote_path{""};
        int local_version{-1};
        int remote_version{-1};
    };

    const std::string default_database_repo_host{"https://api.github.com"};
    const std::vector<std::string> default_database_repo_paths{
        "/repos/pairinteraction/database-sqdt/releases/latest"};

    std::filesystem::path databasedir;
    bool auto_update;
    std::unique_ptr<duckdb::DuckDB> db;
    std::unique_ptr<duckdb::Connection> con;
    std::vector<httplib::Client> pool;
    std::unordered_map<std::string, Table> tables;

    void ensure_presence_of_table(std::string name);
    void ensure_quantum_number_n_is_allowed(std::string name);
};

extern template KetAtom<float>
Database::get_ket<float>(std::string species,
                         const AtomDescriptionByParameters<float> &description);
extern template KetAtom<double>
Database::get_ket<double>(std::string species,
                          const AtomDescriptionByParameters<double> &description);
extern template BasisAtom<float>
Database::get_basis<float>(std::string species, const AtomDescriptionByRanges<float> &description,
                           std::vector<size_t> additional_ket_ids);
extern template BasisAtom<double>
Database::get_basis<double>(std::string species, const AtomDescriptionByRanges<double> &description,
                            std::vector<size_t> additional_ket_ids);
extern template BasisAtom<std::complex<float>> Database::get_basis<std::complex<float>>(
    std::string species, const AtomDescriptionByRanges<std::complex<float>> &description,
    std::vector<size_t> additional_ket_ids);
extern template BasisAtom<std::complex<double>> Database::get_basis<std::complex<double>>(
    std::string species, const AtomDescriptionByRanges<std::complex<double>> &description,
    std::vector<size_t> additional_ket_ids);
extern template Eigen::SparseMatrix<float, Eigen::RowMajor>
Database::get_operator<float>(const BasisAtom<float> &basis, OperatorType type, int q);
extern template Eigen::SparseMatrix<double, Eigen::RowMajor>
Database::get_operator<double>(const BasisAtom<double> &basis, OperatorType type, int q);
extern template Eigen::SparseMatrix<std::complex<float>, Eigen::RowMajor>
Database::get_operator<std::complex<float>>(const BasisAtom<std::complex<float>> &basis,
                                            OperatorType type, int q);
extern template Eigen::SparseMatrix<std::complex<double>, Eigen::RowMajor>
Database::get_operator<std::complex<double>>(const BasisAtom<std::complex<double>> &basis,
                                             OperatorType type, int q);
