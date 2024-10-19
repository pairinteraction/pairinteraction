#include "pairinteraction/basis/BasisCombined.hpp"

#include "pairinteraction/basis/BasisAtom.hpp"
#include "pairinteraction/ket/KetAtom.hpp"
#include "pairinteraction/ket/KetCombined.hpp"
#include "pairinteraction/utils/Range.hpp"

#include <memory>
#include <vector>

namespace pairinteraction {
template <typename Scalar>
BasisCombined<Scalar>::BasisCombined(Private /*unused*/, ketvec_t &&kets, std::string &&id_of_kets,
                                     map_range_t &&map_range_of_state_index2,
                                     std::shared_ptr<const BasisAtom<Scalar>> basis1,
                                     std::shared_ptr<const BasisAtom<Scalar>> basis2)
    : Basis<BasisCombined<Scalar>>(std::move(kets), std::move(id_of_kets)),
      map_range_of_state_index2(std::move(map_range_of_state_index2)), basis1(std::move(basis1)),
      basis2(std::move(basis2)) {}

template <typename Scalar>
const typename BasisCombined<Scalar>::range_t &
BasisCombined<Scalar>::get_index_range(size_t state_index1) const {
    return map_range_of_state_index2.at(state_index1);
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
