#include "pairinteraction/basis/BasisCombined.hpp"

#include "pairinteraction/basis/BasisAtom.hpp"
#include "pairinteraction/ket/KetAtom.hpp"
#include "pairinteraction/ket/KetCombined.hpp"
#include "pairinteraction/utils/Range.hpp"

#include <memory>
#include <vector>

namespace pairinteraction {
template <typename Scalar>
BasisCombined<Scalar>::BasisCombined(Private /*unused*/, ketvec_t &&kets,
                                     map_size_t &&map_index_combined_state,
                                     map_range_t &&map_range_of_index_state2,
                                     size_t number_of_index_state2)
    : Basis<BasisCombined<Scalar>>(std::move(kets)),
      map_index_combined_state(std::move(map_index_combined_state)),
      map_range_of_index_state2(std::move(map_range_of_index_state2)),
      number_of_index_state2(number_of_index_state2) {}

template <typename Scalar>
const typename BasisCombined<Scalar>::range_t &
BasisCombined<Scalar>::get_index_range(size_t index_state1) const {
    return map_range_of_index_state2.at(index_state1);
}

template <typename Scalar>
size_t BasisCombined<Scalar>::get_combined_index(size_t index_state1, size_t index_state2) const {
    return map_index_combined_state.at(index_state1 * number_of_index_state2 + index_state2);
}

// Explicit instantiations
template class BasisCombined<float>;
template class BasisCombined<double>;
template class BasisCombined<std::complex<float>>;
template class BasisCombined<std::complex<double>>;
} // namespace pairinteraction
