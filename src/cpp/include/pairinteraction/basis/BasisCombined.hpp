#pragma once

#include "pairinteraction/basis/Basis.hpp"
#include "pairinteraction/utils/Range.hpp"
#include "pairinteraction/utils/traits.hpp"

#include <complex>
#include <unordered_map>

namespace pairinteraction {
class Database;

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

    using Type = BasisCombined<Scalar>;
    using real_t = typename traits::CrtpTraits<Type>::real_t;
    using ketvec_t = typename traits::CrtpTraits<Type>::ketvec_t;

public:
    using range_t = Range<size_t>;
    using map_size_t = std::unordered_map<size_t, size_t>;
    using map_range_t = std::unordered_map<size_t, range_t>;

    BasisCombined(ketvec_t &&kets, map_size_t &&map_index_combined_state,
                  map_range_t &&map_index_range, size_t number_of_index_state2);
    const range_t &get_index_range(size_t index_state1) const;
    size_t get_combined_index(size_t index_state1, size_t index_state2) const;

private:
    map_size_t map_index_combined_state;
    map_range_t map_range_of_index_state2;
    size_t number_of_index_state2;
};

extern template class BasisCombined<float>;
extern template class BasisCombined<double>;
extern template class BasisCombined<std::complex<float>>;
extern template class BasisCombined<std::complex<double>>;
} // namespace pairinteraction

// TODO generalize this class to an arbitrary number of constituents, also for KetCombined (it
// should take a vector<Ket> for the kets with max overlap)
// TODO create a factor instead having to pass map_index_combined_state etc. to the constructor
// TODO make the constructor private via the Private struct idiom, also for KetCombined
