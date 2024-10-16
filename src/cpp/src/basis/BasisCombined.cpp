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
                                     std::shared_ptr<const BasisAtom<Scalar>> basis1,
                                     std::shared_ptr<const BasisAtom<Scalar>> basis2)
    : Basis<BasisCombined<Scalar>>(std::move(kets)),
      map_index_combined_state(std::move(map_index_combined_state)),
      map_range_of_index_state2(std::move(map_range_of_index_state2)), basis1(std::move(basis1)),
      basis2(std::move(basis2)) {}

template <typename Scalar>
bool BasisCombined<Scalar>::are_valid_indices(size_t index_state1, size_t index_state2) const {
    return map_index_combined_state.count(index_state1 * basis2->get_number_of_states() +
                                          index_state2) > 0;
}

template <typename Scalar>
const typename BasisCombined<Scalar>::range_t &
BasisCombined<Scalar>::get_index_range(size_t index_state1) const {
    return map_range_of_index_state2.at(index_state1);
}

template <typename Scalar>
size_t BasisCombined<Scalar>::get_combined_index(size_t index_state1, size_t index_state2) const {
    return map_index_combined_state.at(index_state1 * basis2->get_number_of_states() +
                                       index_state2);
}

template <typename Scalar>
std::shared_ptr<const BasisAtom<Scalar>> BasisCombined<Scalar>::get_basis1() const {
    return basis1;
}

template <typename Scalar>
std::shared_ptr<const BasisAtom<Scalar>> BasisCombined<Scalar>::get_basis2() const {
    return basis2;
}

// Explicit instantiations
template class BasisCombined<float>;
template class BasisCombined<double>;
template class BasisCombined<std::complex<float>>;
template class BasisCombined<std::complex<double>>;
} // namespace pairinteraction
