#pragma once

#include <filesystem>
#include <memory>
#include <optional>
#include <string>
#include <unordered_map>
#include <vector>

template <typename Real>
class KetAtom;

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

    template <typename Real>
    struct KetsResult {
        std::vector<std::shared_ptr<const KetAtom<Real>>> kets;
        std::string table;
    };

    Database(bool auto_update = false);
    ~Database();
    std::vector<AvailabilitySpecies> get_availability_of_species();
    AvailabilityWigner get_availability_of_wigner_table();

    template <typename Real>
    KetAtom<Real>
    get_ket(std::string species, std::optional<Real> energy, std::optional<float> quantum_number_f,
            std::optional<float> quantum_number_m, std::optional<int> parity,
            std::optional<int> quantum_number_n, std::optional<Real> quantum_number_nu,
            std::optional<Real> quantum_number_l, std::optional<Real> quantum_number_s,
            std::optional<Real> quantum_number_j);

    template <typename Real>
    KetsResult<Real>
    get_kets(std::string species, std::optional<Real> min_energy, std::optional<Real> max_energy,
             std::optional<float> min_quantum_number_f, std::optional<float> max_quantum_number_f,
             std::optional<float> min_quantum_number_m, std::optional<float> max_quantum_number_m,
             std::optional<int> parity, std::optional<int> min_quantum_number_n,
             std::optional<int> max_quantum_number_n, std::optional<Real> min_quantum_number_nu,
             std::optional<Real> max_quantum_number_nu, std::optional<Real> min_quantum_number_l,
             std::optional<Real> max_quantum_number_l, std::optional<Real> min_quantum_number_s,
             std::optional<Real> max_quantum_number_s, std::optional<Real> min_quantum_number_j,
             std::optional<Real> max_quantum_number_j, std::vector<size_t> additional_ket_ids);

private:
    struct Table {
        std::filesystem::path local_path{""};
        std::filesystem::path remote_path{""};
        int local_version{-1};
        int remote_version{-1};
    };

    const std::vector<std::string> database_repo_endpoints{
        "/repos/pairinteraction/database-mqdt/releases/latest",
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

extern template KetAtom<float> Database::get_ket<float>(
    std::string species, std::optional<float> energy, std::optional<float> quantum_number_f,
    std::optional<float> quantum_number_m, std::optional<int> parity,
    std::optional<int> quantum_number_n, std::optional<float> quantum_number_nu,
    std::optional<float> quantum_number_l, std::optional<float> quantum_number_s,
    std::optional<float> quantum_number_j);
extern template KetAtom<double> Database::get_ket<double>(
    std::string species, std::optional<double> energy, std::optional<float> quantum_number_f,
    std::optional<float> quantum_number_m, std::optional<int> parity,
    std::optional<int> quantum_number_n, std::optional<double> quantum_number_nu,
    std::optional<double> quantum_number_l, std::optional<double> quantum_number_s,
    std::optional<double> quantum_number_j);
extern template Database::KetsResult<float> Database::get_kets<float>(
    std::string species, std::optional<float> min_energy, std::optional<float> max_energy,
    std::optional<float> min_quantum_number_f, std::optional<float> max_quantum_number_f,
    std::optional<float> min_quantum_number_m, std::optional<float> max_quantum_number_m,
    std::optional<int> parity, std::optional<int> min_quantum_number_n,
    std::optional<int> max_quantum_number_n, std::optional<float> min_quantum_number_nu,
    std::optional<float> max_quantum_number_nu, std::optional<float> min_quantum_number_l,
    std::optional<float> max_quantum_number_l, std::optional<float> min_quantum_number_s,
    std::optional<float> max_quantum_number_s, std::optional<float> min_quantum_number_j,
    std::optional<float> max_quantum_number_j, std::vector<size_t> additional_ket_ids);
extern template Database::KetsResult<double> Database::get_kets<double>(
    std::string species, std::optional<double> min_energy, std::optional<double> max_energy,
    std::optional<float> min_quantum_number_f, std::optional<float> max_quantum_number_f,
    std::optional<float> min_quantum_number_m, std::optional<float> max_quantum_number_m,
    std::optional<int> parity, std::optional<int> min_quantum_number_n,
    std::optional<int> max_quantum_number_n, std::optional<double> min_quantum_number_nu,
    std::optional<double> max_quantum_number_nu, std::optional<double> min_quantum_number_l,
    std::optional<double> max_quantum_number_l, std::optional<double> min_quantum_number_s,
    std::optional<double> max_quantum_number_s, std::optional<double> min_quantum_number_j,
    std::optional<double> max_quantum_number_j, std::vector<size_t> additional_ket_ids);
