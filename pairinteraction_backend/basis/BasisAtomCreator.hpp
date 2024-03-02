#pragma once

#include "basis/BasisAtom.hpp"
#include "ket/KetAtom.hpp"
#include "utils/Traits.hpp"
#include <complex>
#include <optional>
#include <string>

/**
 * @class BasisAtomCreator
 *
 * @brief Builder class for creating BasisAtom objects.
 *
 * @tparam Scalar Complex number type.
 */
template <typename Scalar>
class BasisAtomCreator {
public:
    using real_t = typename internal::NumTraits<Scalar>::real_t;
    using ket_t = KetAtom<real_t>;
    BasisAtomCreator(std::string species);
    BasisAtomCreator<Scalar> &restrict_energy(Scalar min, Scalar max);
    BasisAtomCreator<Scalar> &restrict_quantum_number_f(float min, float max);
    BasisAtomCreator<Scalar> &restrict_quantum_number_m(float min, float max);
    BasisAtomCreator<Scalar> &restrict_parity(int parity);
    BasisAtomCreator<Scalar> &restrict_quantum_number_n(int min, int max);
    BasisAtomCreator<Scalar> &restrict_quantum_number_nu(Scalar min, Scalar max);
    BasisAtomCreator<Scalar> &restrict_quantum_number_l(Scalar min, Scalar max);
    BasisAtomCreator<Scalar> &restrict_quantum_number_s(Scalar min, Scalar max);
    BasisAtomCreator<Scalar> &restrict_quantum_number_j(Scalar min, Scalar max);
    BasisAtom<Scalar> create() const;

private:
    std::string species;
    std::optional<Scalar> min_energy;
    std::optional<Scalar> max_energy;
    std::optional<float> min_quantum_number_f;
    std::optional<float> max_quantum_number_f;
    std::optional<float> min_quantum_number_m;
    std::optional<float> max_quantum_number_m;
    std::optional<int> parity;
    std::optional<int> min_quantum_number_n;
    std::optional<int> max_quantum_number_n;
    std::optional<Scalar> min_quantum_number_nu;
    std::optional<Scalar> max_quantum_number_nu;
    std::optional<Scalar> min_quantum_number_l;
    std::optional<Scalar> max_quantum_number_l;
    std::optional<Scalar> min_quantum_number_s;
    std::optional<Scalar> max_quantum_number_s;
    std::optional<Scalar> min_quantum_number_j;
    std::optional<Scalar> max_quantum_number_j;
};

extern template class BasisAtomCreator<float>;
extern template class BasisAtomCreator<double>;
extern template class BasisAtomCreator<std::complex<float>>;
extern template class BasisAtomCreator<std::complex<double>>;
