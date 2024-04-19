#include "database/Database.hpp"
#include "ket/KetAtom.hpp"
#include "ket/KetAtomCreator.hpp"

#include <duckdb.hpp>
#include <httplib.h>
#include <nlohmann/json.hpp>

template <typename Real>
KetAtom<Real>
Database::get_ket(std::string species, std::optional<Real> energy,
                  std::optional<float> quantum_number_f, std::optional<float> quantum_number_m,
                  std::optional<int> parity, std::optional<int> quantum_number_n,
                  std::optional<Real> quantum_number_nu, std::optional<Real> quantum_number_l,
                  std::optional<Real> quantum_number_s, std::optional<Real> quantum_number_j) {

    // TODO perform database request and delete the following mocked code

    auto ket =
        KetAtom<Real>(energy.value_or(std::numeric_limits<Real>::quiet_NaN()),
                      quantum_number_f.value_or(std::numeric_limits<float>::quiet_NaN()),
                      quantum_number_m.value_or(std::numeric_limits<float>::quiet_NaN()),
                      parity.value_or(-1), "", 1000, species, quantum_number_n.value_or(0),
                      quantum_number_nu.value_or(std::numeric_limits<Real>::quiet_NaN()), 0,
                      quantum_number_l.value_or(std::numeric_limits<Real>::quiet_NaN()), 0,
                      quantum_number_s.value_or(std::numeric_limits<Real>::quiet_NaN()), 0,
                      quantum_number_j.value_or(std::numeric_limits<Real>::quiet_NaN()), 0, *this);

    return ket;
}

template <typename Real>
std::vector<std::shared_ptr<const KetAtom<Real>>> Database::get_kets(
    std::string species, std::optional<Real> min_energy, std::optional<Real> max_energy,
    std::optional<float> min_quantum_number_f, std::optional<float> max_quantum_number_f,
    std::optional<float> min_quantum_number_m, std::optional<float> max_quantum_number_m,
    std::optional<int> parity, std::optional<int> min_quantum_number_n,
    std::optional<int> max_quantum_number_n, std::optional<Real> min_quantum_number_nu,
    std::optional<Real> max_quantum_number_nu, std::optional<Real> min_quantum_number_l,
    std::optional<Real> max_quantum_number_l, std::optional<Real> min_quantum_number_s,
    std::optional<Real> max_quantum_number_s, std::optional<Real> min_quantum_number_j,
    std::optional<Real> max_quantum_number_j) {

    // TODO perform database request and delete the following mocked code

    Real energy = 0;
    float quantum_number_f = 0;
    float quantum_number_m = 0;
    int p = -1;
    size_t id = 1000;
    int quantum_number_n = 0;
    Real quantum_number_nu = 0;
    Real quantum_number_l = 0;
    Real quantum_number_s = 0;
    Real quantum_number_j = 0;

    std::vector<std::shared_ptr<const KetAtom<Real>>> kets;
    kets.push_back(std::make_shared<const KetAtom<Real>>(
        KetAtom<Real>(energy, quantum_number_f, quantum_number_m, p, species, id, species,
                      quantum_number_n, quantum_number_nu, 0, quantum_number_l, 0, quantum_number_s,
                      0, quantum_number_j, 0, *this)));

    return kets;
}

// Explicit instantiations
template KetAtom<float> Database::get_ket<float>(
    std::string species, std::optional<float> energy, std::optional<float> quantum_number_f,
    std::optional<float> quantum_number_m, std::optional<int> parity,
    std::optional<int> quantum_number_n, std::optional<float> quantum_number_nu,
    std::optional<float> quantum_number_l, std::optional<float> quantum_number_s,
    std::optional<float> quantum_number_j);
template KetAtom<double> Database::get_ket<double>(
    std::string species, std::optional<double> energy, std::optional<float> quantum_number_f,
    std::optional<float> quantum_number_m, std::optional<int> parity,
    std::optional<int> quantum_number_n, std::optional<double> quantum_number_nu,
    std::optional<double> quantum_number_l, std::optional<double> quantum_number_s,
    std::optional<double> quantum_number_j);
template std::vector<std::shared_ptr<const KetAtom<float>>> Database::get_kets<float>(
    std::string species, std::optional<float> min_energy, std::optional<float> max_energy,
    std::optional<float> min_quantum_number_f, std::optional<float> max_quantum_number_f,
    std::optional<float> min_quantum_number_m, std::optional<float> max_quantum_number_m,
    std::optional<int> parity, std::optional<int> min_quantum_number_n,
    std::optional<int> max_quantum_number_n, std::optional<float> min_quantum_number_nu,
    std::optional<float> max_quantum_number_nu, std::optional<float> min_quantum_number_l,
    std::optional<float> max_quantum_number_l, std::optional<float> min_quantum_number_s,
    std::optional<float> max_quantum_number_s, std::optional<float> min_quantum_number_j,
    std::optional<float> max_quantum_number_j);
template std::vector<std::shared_ptr<const KetAtom<double>>> Database::get_kets<double>(
    std::string species, std::optional<double> min_energy, std::optional<double> max_energy,
    std::optional<float> min_quantum_number_f, std::optional<float> max_quantum_number_f,
    std::optional<float> min_quantum_number_m, std::optional<float> max_quantum_number_m,
    std::optional<int> parity, std::optional<int> min_quantum_number_n,
    std::optional<int> max_quantum_number_n, std::optional<double> min_quantum_number_nu,
    std::optional<double> max_quantum_number_nu, std::optional<double> min_quantum_number_l,
    std::optional<double> max_quantum_number_l, std::optional<double> min_quantum_number_s,
    std::optional<double> max_quantum_number_s, std::optional<double> min_quantum_number_j,
    std::optional<double> max_quantum_number_j);
