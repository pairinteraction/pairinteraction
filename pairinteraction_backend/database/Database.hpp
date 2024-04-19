#pragma once

#include <memory>
#include <optional>
#include <string>
#include <vector>

template <typename Real>
class KetAtom;

class Database {
public:
    Database() = default;

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
};

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
