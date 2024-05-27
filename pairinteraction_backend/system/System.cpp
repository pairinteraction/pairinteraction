#include "system/System.hpp"
#include "operator/Operator.hpp"
#include "operator/OperatorAtom.hpp"

#include <Eigen/SparseCore>
#include <complex>
#include <memory>

template <typename Derived>
System<Derived>::System(std::shared_ptr<const basis_t> basis)
    : hamiltonian(std::make_shared<typename System<Derived>::operator_t>(basis)) {}

template <typename Derived>
const Derived &System<Derived>::derived() const {
    return static_cast<const Derived &>(*this);
}

template <typename Derived>
std::shared_ptr<const typename System<Derived>::basis_t> System<Derived>::get_basis() const {
    return hamiltonian->get_basis();
}

template <typename Derived>
const Eigen::SparseMatrix<typename System<Derived>::scalar_t, Eigen::RowMajor> &
System<Derived>::get_matrix() const {
    return hamiltonian->get_matrix();
}

template <typename Derived>
const Transformation<typename System<Derived>::scalar_t> &
System<Derived>::get_transformation() const {
    return hamiltonian->get_transformation();
}

template <typename Derived>
Transformation<typename System<Derived>::scalar_t>
System<Derived>::get_rotator(real_t alpha, real_t beta, real_t gamma) const {
    return hamiltonian->get_rotator(alpha, beta, gamma);
}

template <typename Derived>
Sorting System<Derived>::get_sorter(TransformationType label) const {
    return hamiltonian->get_sorter(label);
}

template <typename Derived>
Blocks System<Derived>::get_blocks(TransformationType label) const {
    return hamiltonian->get_blocks(label);
}

template <typename Derived>
Derived System<Derived>::transformed(const Transformation<scalar_t> &transformation) const {
    Derived transformed = derived();
    transformed.hamiltonian =
        std::make_shared<operator_t>(hamiltonian->transformed(transformation));
    return transformed;
}

template <typename Derived>
Derived System<Derived>::transformed(const Sorting &transformation) const {
    Derived transformed = derived();
    transformed.hamiltonian =
        std::make_shared<operator_t>(hamiltonian->transformed(transformation));
    return transformed;
}

// Explicit instantiation
#include "system/SystemAtom.hpp"

template class System<SystemAtom<float>>;
template class System<SystemAtom<double>>;
template class System<SystemAtom<std::complex<float>>>;
template class System<SystemAtom<std::complex<double>>>;
