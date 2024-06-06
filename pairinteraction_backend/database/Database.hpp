#pragma once

#include "utils/traits.hpp"
#include <complex>
#include <filesystem>
#include <memory>
#include <optional>
#include <string>
#include <unordered_map>
#include <vector>

enum class OperatorType;

template <typename Scalar>
struct AtomDescriptionByParameters;

template <typename Scalar>
struct AtomDescriptionByRanges;

template <typename Real>
class KetAtom;

template <typename Scalar>
class BasisAtom;

template <typename Scalar>
class OperatorAtom;

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

    Database(std::optional<bool> download_missing = std::optional<bool>(),
             std::filesystem::path databasedir = "");
    ~Database();
    std::vector<AvailabilitySpecies> get_availability_of_species();
    AvailabilityWigner get_availability_of_wigner_table();
    static Database &
    get_global_instance(std::optional<bool> download_missing = std::optional<bool>(),
                        std::filesystem::path databasedir = "");

    template <typename Real>
    std::shared_ptr<const KetAtom<Real>>
    get_ket(std::string species, const AtomDescriptionByParameters<Real> &description);

    template <typename Scalar>
    std::shared_ptr<const BasisAtom<Scalar>> get_basis(
        std::string species,
        const AtomDescriptionByRanges<typename traits::NumTraits<Scalar>::real_t> &description,
        std::vector<size_t> additional_ket_ids);

    template <typename Scalar>
    OperatorAtom<Scalar> get_operator(std::shared_ptr<const BasisAtom<Scalar>> basis,
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
        "/repos/pairinteraction/database-sqdt/releases/latest",
        "/repos/pairinteraction/database-mqdt/releases/latest"};

    std::filesystem::path databasedir;
    bool download_missing;
    std::unique_ptr<duckdb::DuckDB> db;
    std::unique_ptr<duckdb::Connection> con;
    std::vector<httplib::Client> pool;
    std::unordered_map<std::string, Table> tables;

    void ensure_presence_of_table(std::string name);
    void ensure_quantum_number_n_is_allowed(std::string name);
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
extern template OperatorAtom<float>
Database::get_operator<float>(std::shared_ptr<const BasisAtom<float>> basis, OperatorType type,
                              int q);
extern template OperatorAtom<double>
Database::get_operator<double>(std::shared_ptr<const BasisAtom<double>> basis, OperatorType type,
                               int q);
extern template OperatorAtom<std::complex<float>> Database::get_operator<std::complex<float>>(
    std::shared_ptr<const BasisAtom<std::complex<float>>> basis, OperatorType type, int q);
extern template OperatorAtom<std::complex<double>> Database::get_operator<std::complex<double>>(
    std::shared_ptr<const BasisAtom<std::complex<double>>> basis, OperatorType type, int q);
