#pragma once

#include "pairinteraction/utils/Range.hpp"
#include "pairinteraction/utils/traits.hpp"

#include <complex>
#include <memory>
#include <vector>

namespace pairinteraction {
template <typename Scalar>
class BasisCombined;

template <typename Scalar>
class SystemAtom;

template <typename Scalar>
class KetCombined;

template <typename Scalar>
class BasisCombinedCreator {
    static_assert(traits::NumTraits<Scalar>::from_floating_point_v);

public:
    using real_t = typename traits::NumTraits<Scalar>::real_t;
    using basis_t = BasisCombined<Scalar>;
    using ket_t = KetCombined<real_t>;
    using ketvec_t = std::vector<std::shared_ptr<const ket_t>>;

    BasisCombinedCreator() = default;
    BasisCombinedCreator<Scalar> &add(const SystemAtom<Scalar> &system_atom);
    BasisCombinedCreator<Scalar> &restrict_energy(real_t min, real_t max);
    BasisCombinedCreator<Scalar> &restrict_quantum_number_m(real_t min, real_t max);
    std::shared_ptr<const BasisCombined<Scalar>> create() const;

private:
    std::vector<std::reference_wrapper<const SystemAtom<Scalar>>> systems_atom;
    Range<real_t> range_energy;
    Range<real_t> range_quantum_number_m;
};

extern template class BasisCombinedCreator<float>;
extern template class BasisCombinedCreator<double>;
extern template class BasisCombinedCreator<std::complex<float>>;
extern template class BasisCombinedCreator<std::complex<double>>;
} // namespace pairinteraction
