#include "pintr/system/System.hpp"

#include "pintr/operator/OperatorAtom.hpp"
#include "pintr/system/SystemAtom.hpp"
#include "pintr/utils/eigen_assertion.hpp"

#include <Eigen/SparseCore>
#include <complex>
#include <limits>
#include <memory>
#include <oneapi/tbb.h>

namespace pintr {
template <typename Derived>
System<Derived>::System(std::shared_ptr<const basis_t> basis)
    : hamiltonian(std::make_unique<typename System<Derived>::operator_t>(basis)) {}

template <typename Derived>
System<Derived>::System(const System &other)
    : hamiltonian(std::make_unique<typename System<Derived>::operator_t>(*other.hamiltonian)) {}

template <typename Derived>
System<Derived>::~System() = default;

template <typename Derived>
const Derived &System<Derived>::derived() const {
    return static_cast<const Derived &>(*this);
}

template <typename Derived>
std::shared_ptr<const typename System<Derived>::basis_t> System<Derived>::get_basis() const {
    if (hamiltonian_requires_construction) {
        construct_hamiltonian();
        hamiltonian_requires_construction = false;
    }
    return hamiltonian->get_basis();
}

template <typename Derived>
const Eigen::SparseMatrix<typename System<Derived>::scalar_t, Eigen::RowMajor> &
System<Derived>::get_matrix() const {
    if (hamiltonian_requires_construction) {
        construct_hamiltonian();
        hamiltonian_requires_construction = false;
    }
    return hamiltonian->get_matrix();
}

template <typename Derived>
const Transformation<typename System<Derived>::scalar_t> &
System<Derived>::get_transformation() const {
    if (hamiltonian_requires_construction) {
        construct_hamiltonian();
        hamiltonian_requires_construction = false;
    }
    return hamiltonian->get_transformation();
}

template <typename Derived>
Transformation<typename System<Derived>::scalar_t>
System<Derived>::get_rotator(real_t alpha, real_t beta, real_t gamma) const {
    if (hamiltonian_requires_construction) {
        construct_hamiltonian();
        hamiltonian_requires_construction = false;
    }
    return hamiltonian->get_rotator(alpha, beta, gamma);
}

template <typename Derived>
Sorting System<Derived>::get_sorter(TransformationType label) const {
    if (hamiltonian_requires_construction) {
        construct_hamiltonian();
        hamiltonian_requires_construction = false;
    }
    return hamiltonian->get_sorter(label);
}

template <typename Derived>
Blocks System<Derived>::get_blocks(TransformationType label) const {
    if (hamiltonian_requires_construction) {
        construct_hamiltonian();
        hamiltonian_requires_construction = false;
    }
    return hamiltonian->get_blocks(label);
}

template <typename Derived>
Derived System<Derived>::transformed(const Transformation<scalar_t> &transformation) const {
    if (hamiltonian_requires_construction) {
        construct_hamiltonian();
        hamiltonian_requires_construction = false;
    }
    Derived transformed = derived();
    transformed.hamiltonian =
        std::make_unique<operator_t>(hamiltonian->transformed(transformation));
    return transformed;
}

template <typename Derived>
Derived System<Derived>::transformed(const Sorting &transformation) const {
    if (hamiltonian_requires_construction) {
        construct_hamiltonian();
        hamiltonian_requires_construction = false;
    }
    Derived transformed = derived();
    transformed.hamiltonian =
        std::make_unique<operator_t>(hamiltonian->transformed(transformation));
    return transformed;
}

template <typename Derived>
void System<Derived>::diagonalize(const DiagonalizerInterface<scalar_t> &diagonalizer,
                                  int precision) const {
    real_t min_eigenvalue = std::numeric_limits<real_t>::min();
    real_t max_eigenvalue = std::numeric_limits<real_t>::max();
    diagonalize(diagonalizer, min_eigenvalue, max_eigenvalue, precision);
}

template <typename Derived>
void System<Derived>::diagonalize(const DiagonalizerInterface<scalar_t> &diagonalizer,
                                  real_t min_eigenvalue, real_t max_eigenvalue,
                                  int precision) const {
    if (hamiltonian_requires_construction) {
        construct_hamiltonian();
        hamiltonian_requires_construction = false;
    }

    // TODO get TransformationType from derived class and diagonalize blocks in parallel
    bool Data[1000];
    oneapi::tbb::parallel_for(0, 1000, 1, [&Data](int i) { Data[i] = true; });

    auto eigensys =
        diagonalizer.eigh(hamiltonian->get_matrix(), min_eigenvalue, max_eigenvalue, precision);

    // Store the diagonalized hamiltonian (possible future optimization: use
    // eigensys.eigenvalues directly instead of transforming the hamiltonian with the eigenvectors)
    hamiltonian = std::make_unique<operator_t>(hamiltonian->transformed(eigensys.eigenvectors));
}

// Explicit instantiation
template class System<SystemAtom<float>>;
template class System<SystemAtom<double>>;
template class System<SystemAtom<std::complex<float>>>;
template class System<SystemAtom<std::complex<double>>>;
} // namespace pintr
