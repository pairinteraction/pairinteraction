#include "pairinteraction/operator/OperatorPair.hpp"

#include "pairinteraction/basis/BasisPair.hpp"
#include "pairinteraction/enums/OperatorType.hpp"

namespace pairinteraction {

template <typename Scalar>
OperatorPair<Scalar>::OperatorPair(std::shared_ptr<const basis_t> basis)
    : Operator<OperatorPair<Scalar>>(std::move(basis)) {}

template <typename Scalar>
OperatorPair<Scalar>::OperatorPair(std::shared_ptr<const basis_t> basis, OperatorType type)
    : Operator<OperatorPair<Scalar>>(std::move(basis)) {
    if (type == OperatorType::ENERGY) {
        this->initialize_as_energy_operator();
    } else {
        throw std::invalid_argument("Only OperatorType::ENERGY is supported.");
    }
}

// Explicit instantiations
template class OperatorPair<float>;
template class OperatorPair<double>;
template class OperatorPair<std::complex<float>>;
template class OperatorPair<std::complex<double>>;
} // namespace pairinteraction
