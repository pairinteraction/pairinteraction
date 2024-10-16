#include "pairinteraction/operator/OperatorCombined.hpp"

#include "pairinteraction/basis/BasisCombined.hpp"
#include "pairinteraction/enums/OperatorType.hpp"

namespace pairinteraction {

template <typename Scalar>
OperatorCombined<Scalar>::OperatorCombined(std::shared_ptr<const basis_t> basis)
    : Operator<OperatorCombined<Scalar>>(std::move(basis)) {}

template <typename Scalar>
OperatorCombined<Scalar>::OperatorCombined(std::shared_ptr<const basis_t> basis, OperatorType type)
    : Operator<OperatorCombined<Scalar>>(std::move(basis)) {
    if (type == OperatorType::ENERGY) {
        this->initialize_as_energy_operator();
    } else {
        throw std::invalid_argument("Only OperatorType::ENERGY is supported.");
    }
}

// Explicit instantiations
template class OperatorCombined<float>;
template class OperatorCombined<double>;
template class OperatorCombined<std::complex<float>>;
template class OperatorCombined<std::complex<double>>;
} // namespace pairinteraction
