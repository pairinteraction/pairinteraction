#pragma once

#include "pintr/enums/Parity.hpp"
#include "pintr/utils/Range.hpp"
#include "pintr/utils/traits.hpp"

#include <complex>
#include <memory>
#include <optional>
#include <string>
#include <vector>

template <typename Scalar>
class BasisAtom;

template <typename Real>
class KetAtom;

class Database;

/**
 * @class BasisAtomCreator
 *
 * @brief Builder class for creating BasisAtom objects.
 *
 * @tparam Scalar Complex number type.
 */
template <typename Scalar>
class BasisAtomCreator {
    static_assert(traits::NumTraits<Scalar>::from_floating_point_v);

public:
    using real_t = typename traits::NumTraits<Scalar>::real_t;
    using ket_t = KetAtom<real_t>;
    BasisAtomCreator() = default;
    BasisAtomCreator<Scalar> &set_species(const std::string &value);
    BasisAtomCreator<Scalar> &restrict_energy(real_t min, real_t max);
    BasisAtomCreator<Scalar> &restrict_quantum_number_f(real_t min, real_t max);
    BasisAtomCreator<Scalar> &restrict_quantum_number_m(real_t min, real_t max);
    BasisAtomCreator<Scalar> &restrict_parity(Parity value);
    BasisAtomCreator<Scalar> &restrict_quantum_number_n(int min, int max);
    BasisAtomCreator<Scalar> &restrict_quantum_number_nu(real_t min, real_t max);
    BasisAtomCreator<Scalar> &restrict_quantum_number_l(real_t min, real_t max);
    BasisAtomCreator<Scalar> &restrict_quantum_number_s(real_t min, real_t max);
    BasisAtomCreator<Scalar> &restrict_quantum_number_j(real_t min, real_t max);
    BasisAtomCreator<Scalar> &add_ket(std::shared_ptr<const ket_t> ket);
    std::shared_ptr<const BasisAtom<Scalar>> create(Database &database) const;

private:
    std::optional<std::string> species;
    Parity parity{Parity::UNKNOWN};
    Range<real_t> range_energy;
    Range<real_t> range_quantum_number_f;
    Range<real_t> range_quantum_number_m;
    Range<int> range_quantum_number_n;
    Range<real_t> range_quantum_number_nu;
    Range<real_t> range_quantum_number_l;
    Range<real_t> range_quantum_number_s;
    Range<real_t> range_quantum_number_j;
    std::vector<size_t> additional_ket_ids;
    std::optional<std::string> additional_ket_species;
};

extern template class BasisAtomCreator<float>;
extern template class BasisAtomCreator<double>;
extern template class BasisAtomCreator<std::complex<float>>;
extern template class BasisAtomCreator<std::complex<double>>;
