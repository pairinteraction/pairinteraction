#pragma once

#include "enums/OperatorType.hpp"
#include "operator/Operator.hpp"
#include "utils/traits.hpp"

#include <memory>
#include <vector>

class Database;

template <typename Scalar>
class BasisAtom;

template <typename Real>
class KetAtom;

// Specialize OperatorTraits for OperatorAtom
template <typename T>
class OperatorAtom;

template <typename Scalar>
struct traits::OperatorTraits<OperatorAtom<Scalar>> {
    using scalar_t = Scalar;
    using real_t = typename traits::NumTraits<Scalar>::real_t;
    using ket_t = KetAtom<real_t>;
    using ketvec_t = std::vector<std::shared_ptr<const ket_t>>;
    using basis_t = BasisAtom<scalar_t>;
};

template <typename Scalar>
class OperatorAtom : public Operator<OperatorAtom<Scalar>> {
public:
    using Type = OperatorAtom<Scalar>;
    using basis_t = typename traits::OperatorTraits<Type>::basis_t;

    OperatorAtom(const basis_t &basis, OperatorType type, int q);

private:
    friend class Database;
    OperatorType type;
    int q;
};

extern template class OperatorAtom<float>;
extern template class OperatorAtom<double>;
extern template class OperatorAtom<std::complex<float>>;
extern template class OperatorAtom<std::complex<double>>;
