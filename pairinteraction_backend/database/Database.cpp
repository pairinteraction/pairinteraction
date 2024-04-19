#include "database/Database.hpp"
#include "ket/KetAtom.hpp"
#include "ket/KetAtomCreator.hpp"

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

    // TODO perform database request

    std::vector<std::shared_ptr<const KetAtom<Real>>> kets;
    kets.reserve(2);
    kets.push_back(std::make_shared<const KetAtom<Real>>(
        KetAtomCreator<Real>(species, 60, 1, 0.5, -0.5).create(*this)));
    kets.push_back(std::make_shared<const KetAtom<Real>>(
        KetAtomCreator<Real>(species, 60, 1, 0.5, 0.5).create(*this)));

    return kets;
}

// Explicit instantiations
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
