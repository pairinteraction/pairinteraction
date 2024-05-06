#pragma once

#include <limits>
#include <memory>
#include <optional>
#include <string>

class Database;

template <typename Real>
class KetAtom;

/**
 * @class KetAtomCreator
 *
 * @brief Builder class for creating KetAtom objects.
 *
 * @tparam Real Real number type.
 */
template <typename Real>
class KetAtomCreator {
public:
    KetAtomCreator() = default;
    KetAtomCreator(std::string species, int n, Real l, float j, float m);
    KetAtomCreator<Real> &set_species(std::string value);
    KetAtomCreator<Real> &set_energy(Real value);
    KetAtomCreator<Real> &set_quantum_number_f(float value);
    KetAtomCreator<Real> &set_quantum_number_m(float value);
    KetAtomCreator<Real> &set_parity(int value);
    KetAtomCreator<Real> &set_quantum_number_n(int value);
    KetAtomCreator<Real> &set_quantum_number_nu(Real value);
    KetAtomCreator<Real> &set_quantum_number_l(Real value);
    KetAtomCreator<Real> &set_quantum_number_s(Real value);
    KetAtomCreator<Real> &set_quantum_number_j(Real value);
    std::shared_ptr<const KetAtom<Real>> create(Database &database) const;

private:
    std::optional<std::string> species;
    std::optional<Real> energy;
    std::optional<float> quantum_number_f;
    std::optional<float> quantum_number_m;
    std::optional<int> parity;
    std::optional<int> quantum_number_n;
    std::optional<Real> quantum_number_nu;
    std::optional<Real> quantum_number_l;
    std::optional<Real> quantum_number_s;
    std::optional<Real> quantum_number_j;
};

extern template class KetAtomCreator<float>;
extern template class KetAtomCreator<double>;
