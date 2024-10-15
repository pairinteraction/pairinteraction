#pragma once

#include "pairinteraction/basis/Basis.hpp"
#include "pairinteraction/utils/Range.hpp"
#include "pairinteraction/utils/traits.hpp"

#include <complex>
#include <unordered_map>

namespace pairinteraction {
template <typename Scalar>
class SystemCombined;

template <typename Real>
class KetCombined;

template <typename Scalar>
class BasisCombined;

template <typename Scalar>
class BasisAtom;

template <typename Real>
class KetAtom;

template <typename Scalar>
struct traits::CrtpTraits<BasisCombined<Scalar>> {
    using scalar_t = Scalar;
    using real_t = typename traits::NumTraits<Scalar>::real_t;
    using ket_t = KetCombined<real_t>;
    using ketvec_t = std::vector<std::shared_ptr<const ket_t>>;
};

template <typename Scalar>
class BasisCombined : public Basis<BasisCombined<Scalar>> {
    static_assert(traits::NumTraits<Scalar>::from_floating_point_v);

    friend class SystemCombined<Scalar>;
    struct Private {};

public:
    using Type = BasisCombined<Scalar>;
    using real_t = typename traits::CrtpTraits<Type>::real_t;
    using ketvec_t = typename traits::CrtpTraits<Type>::ketvec_t;
    using range_t = Range<size_t>;
    using map_size_t = std::unordered_map<size_t, size_t>;
    using map_range_t = std::unordered_map<size_t, range_t>;

    BasisCombined(Private /*unused*/, ketvec_t &&kets, map_size_t &&map_index_combined_state,
                  map_range_t &&map_index_range, size_t number_of_index_state2);

private:
    const range_t &get_index_range(size_t index_state1) const;
    size_t get_combined_index(size_t index_state1, size_t index_state2) const;

    map_size_t map_index_combined_state;
    map_range_t map_range_of_index_state2;
    size_t number_of_index_state2;
};

extern template class BasisCombined<float>;
extern template class BasisCombined<double>;
extern template class BasisCombined<std::complex<float>>;
extern template class BasisCombined<std::complex<double>>;
} // namespace pairinteraction
