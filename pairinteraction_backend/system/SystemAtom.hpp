#pragma once

#include "system/System.hpp"
#include "utils/traits.hpp"

#include <memory>
#include <vector>

template <typename Scalar>
class OperatorAtom;

template <typename Scalar>
class BasisAtom;

template <typename Real>
class KetAtom;

// Specialize CrtpTraits for SystemAtom
template <typename T>
class SystemAtom;

template <typename Scalar>
struct traits::CrtpTraits<SystemAtom<Scalar>> {
    using scalar_t = Scalar;
    using real_t = typename traits::NumTraits<Scalar>::real_t;
    using ket_t = KetAtom<real_t>;
    using ketvec_t = std::vector<std::shared_ptr<const ket_t>>;
    using basis_t = BasisAtom<scalar_t>;
    using operator_t = OperatorAtom<scalar_t>;
};

template <typename Scalar>
class SystemAtom : public System<SystemAtom<Scalar>> {
public:
    static_assert(traits::NumTraits<Scalar>::from_floating_point_v);

    using Type = SystemAtom<Scalar>;
    using basis_t = typename traits::CrtpTraits<Type>::basis_t;

    SystemAtom(std::shared_ptr<const basis_t> basis);

private:
    void construct_hamiltonian() const override{};
};

extern template class SystemAtom<float>;
extern template class SystemAtom<double>;
extern template class SystemAtom<std::complex<float>>;
extern template class SystemAtom<std::complex<double>>;
