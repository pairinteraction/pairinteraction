#include "operator/Operator.hpp"
#include "basis/BasisAtom.hpp"

#include <numeric>
#include <set>

template <typename Derived>
Operator<Derived>::Operator(const basis_t &basis) : basis(std::make_shared<basis_t>(basis)) {}

template <typename Derived>
const Derived &Operator<Derived>::derived() const {
    return static_cast<const Derived &>(*this);
}

template <typename Derived>
const typename Operator<Derived>::basis_t &Operator<Derived>::get_basis() const {
    return *basis;
}

template <typename Derived>
const Eigen::SparseMatrix<typename Operator<Derived>::scalar_t, Eigen::RowMajor> &
Operator<Derived>::get_matrix() const {
    return matrix;
}

template <typename Derived>
const typename Operator<Derived>::ketvec_t &Operator<Derived>::get_kets() const {
    return basis->get_kets();
}

template <typename Derived>
const Eigen::SparseMatrix<typename Operator<Derived>::scalar_t, Eigen::RowMajor> &
Operator<Derived>::get_coefficients() const {
    return basis->get_coefficients();
}

template <typename Derived>
size_t Operator<Derived>::get_number_of_states() const {
    return basis->get_number_of_states();
}

template <typename Derived>
size_t Operator<Derived>::get_number_of_kets() const {
    return basis->get_number_of_kets();
}

template <typename Derived>
Eigen::SparseMatrix<typename Operator<Derived>::scalar_t>
Operator<Derived>::impl_get_rotator(real_t alpha, real_t beta, real_t gamma) const {
    return basis->get_rotator(alpha, beta, gamma);
}

template <typename Derived>
std::vector<int> Operator<Derived>::impl_get_sorter(SortBy label) const {
    auto sorter = basis->get_sorter(label & ~SortBy::DIAGONAL);

    std::vector<real_t> energies_of_states;
    energies_of_states.reserve(matrix.rows());
    for (int i = 0; i < matrix.rows(); ++i) {
        if constexpr (traits::is_complex_v<scalar_t>) {
            energies_of_states.push_back(matrix.coeff(i, i).real());
        } else {
            energies_of_states.push_back(matrix.coeff(i, i));
        }
    }

    if ((label & SortBy::DIAGONAL) == SortBy::DIAGONAL) {
        std::stable_sort(sorter.begin(), sorter.end(), [&](int i, int j) {
            return energies_of_states[i] < energies_of_states[j];
        });
    }

    return sorter;
}

template <typename Derived>
std::vector<int> Operator<Derived>::impl_get_blocks(SortBy label) const {
    // TODO implement block extraction
    if ((label & SortBy::DIAGONAL) == SortBy::DIAGONAL) {
        throw std::runtime_error("Extracting blocks based on diagonal elements is not supported.");
    }
    return basis->get_blocks(label);
}

template <typename Derived>
void Operator<Derived>::impl_transform(const Eigen::SparseMatrix<scalar_t> &transformator) {
    matrix = transformator.adjoint() * matrix * transformator;
    basis->transform(transformator);
}

template <typename Derived>
void Operator<Derived>::impl_sort(const std::vector<int> &sorter) {
    // TODO sort the matrix
    basis->sort(sorter);
}

// Explicit instantiations
#include "operator/OperatorAtom.hpp"

template class Operator<OperatorAtom<float>>;
template class Operator<OperatorAtom<double>>;
template class Operator<OperatorAtom<std::complex<float>>>;
template class Operator<OperatorAtom<std::complex<double>>>;
