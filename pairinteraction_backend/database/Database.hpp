#pragma once

#include <filesystem>
#include <memory>
#include <optional>
#include <string>
#include <unordered_map>
#include <vector>

template <typename Real>
class KetAtom;

class Database {
public:
    Database(bool auto_update = true);

    template <typename Real>
    KetAtom<Real>
    get_ket(std::string species, std::optional<Real> energy, std::optional<float> quantum_number_f,
            std::optional<float> quantum_number_m, std::optional<int> parity,
            std::optional<int> quantum_number_n, std::optional<Real> quantum_number_nu,
            std::optional<Real> quantum_number_l, std::optional<Real> quantum_number_s,
            std::optional<Real> quantum_number_j);

    template <typename Real>
    std::vector<std::shared_ptr<const KetAtom<Real>>>
    get_kets(std::string species, std::optional<Real> min_energy, std::optional<Real> max_energy,
             std::optional<float> min_quantum_number_f, std::optional<float> max_quantum_number_f,
             std::optional<float> min_quantum_number_m, std::optional<float> max_quantum_number_m,
             std::optional<int> parity, std::optional<int> min_quantum_number_n,
             std::optional<int> max_quantum_number_n, std::optional<Real> min_quantum_number_nu,
             std::optional<Real> max_quantum_number_nu, std::optional<Real> min_quantum_number_l,
             std::optional<Real> max_quantum_number_l, std::optional<Real> min_quantum_number_s,
             std::optional<Real> max_quantum_number_s, std::optional<Real> min_quantum_number_j,
             std::optional<Real> max_quantum_number_j);

private:
    struct LocalDatabase {
        std::filesystem::path path{""};
        int version{-1};
    };
    struct RemoteDatabase {
        std::string url{""};
        int version{-1};
    };
    const std::vector<std::string> database_repo_endpoints{
        {"/repos/MultiQuantum/extending-pairinteraction-proposal/releases/latest",
         "/repos/MultiQuantum/rubidium-pairinteraction-database/releases/latest",
         "/repos/MultiQuantum/strontium-pairinteraction-database/releases/latest"}};
    const std::string github_access_token{
        "github_pat_11ACL5QGA0vJHuTHWEnkex_"
        "hTQiAbwt2R8gLX92pUwR9GofZXq5ee4a1qfhX6SWSEpQQAL2CD3NRQRpm3H"};
    std::filesystem::path databasedir;
    bool auto_update;
    std::unordered_map<std::string, RemoteDatabase> remote_databases;
    std::unordered_map<std::string, LocalDatabase> local_databases;
    LocalDatabase get_local_database(std::string name);
    RemoteDatabase get_remote_database(std::string name);
    void ensure_presence_of_local_database(std::string name);
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
extern template std::vector<std::shared_ptr<const KetAtom<float>>> Database::get_kets<float>(
    std::string species, std::optional<float> min_energy, std::optional<float> max_energy,
    std::optional<float> min_quantum_number_f, std::optional<float> max_quantum_number_f,
    std::optional<float> min_quantum_number_m, std::optional<float> max_quantum_number_m,
    std::optional<int> parity, std::optional<int> min_quantum_number_n,
    std::optional<int> max_quantum_number_n, std::optional<float> min_quantum_number_nu,
    std::optional<float> max_quantum_number_nu, std::optional<float> min_quantum_number_l,
    std::optional<float> max_quantum_number_l, std::optional<float> min_quantum_number_s,
    std::optional<float> max_quantum_number_s, std::optional<float> min_quantum_number_j,
    std::optional<float> max_quantum_number_j);
extern template std::vector<std::shared_ptr<const KetAtom<double>>> Database::get_kets<double>(
    std::string species, std::optional<double> min_energy, std::optional<double> max_energy,
    std::optional<float> min_quantum_number_f, std::optional<float> max_quantum_number_f,
    std::optional<float> min_quantum_number_m, std::optional<float> max_quantum_number_m,
    std::optional<int> parity, std::optional<int> min_quantum_number_n,
    std::optional<int> max_quantum_number_n, std::optional<double> min_quantum_number_nu,
    std::optional<double> max_quantum_number_nu, std::optional<double> min_quantum_number_l,
    std::optional<double> max_quantum_number_l, std::optional<double> min_quantum_number_s,
    std::optional<double> max_quantum_number_s, std::optional<double> min_quantum_number_j,
    std::optional<double> max_quantum_number_j);
