// SPDX-FileCopyrightText: 2024 Pairinteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "pairinteraction/operator/Operator.hpp"

#include "pairinteraction/basis/BasisAtom.hpp"
#include "pairinteraction/basis/BasisPair.hpp"
#include "pairinteraction/enums/OperatorType.hpp"
#include "pairinteraction/enums/TransformationType.hpp"
#include "pairinteraction/ket/KetAtom.hpp"
#include "pairinteraction/ket/KetPair.hpp"
#include "pairinteraction/operator/OperatorAtom.hpp"
#include "pairinteraction/operator/OperatorPair.hpp"
#include "pairinteraction/utils/eigen_assertion.hpp"

#include <Eigen/SparseCore>
#include <memory>

namespace pairinteraction {
template <typename Derived>
Operator<Derived>::Operator(std::shared_ptr<const basis_t> basis) : basis(std::move(basis)) {
    this->matrix = Eigen::SparseMatrix<scalar_t, Eigen::RowMajor>(
        this->basis->get_number_of_states(), this->basis->get_number_of_states());
}

template <typename Derived>
void Operator<Derived>::initialize_as_energy_operator() {
    Eigen::SparseMatrix<scalar_t, Eigen::RowMajor> tmp(this->basis->get_number_of_kets(),
                                                       this->basis->get_number_of_kets());
    tmp.reserve(Eigen::VectorXi::Constant(this->basis->get_number_of_kets(), 1));
    size_t idx = 0;
    for (const auto &ket : this->basis->get_kets()) {
        tmp.insert(idx, idx) = ket->get_energy();
        ++idx;
    }
    tmp.makeCompressed();

    this->matrix =
        this->basis->get_coefficients().adjoint() * tmp * this->basis->get_coefficients();
}

template <typename Derived>
void Operator<Derived>::initialize_from_matrix(
    Eigen::SparseMatrix<scalar_t, Eigen::RowMajor> &&matrix) {
    if (static_cast<size_t>(matrix.rows()) != this->basis->get_number_of_states() ||
        static_cast<size_t>(matrix.cols()) != this->basis->get_number_of_states()) {
        throw std::invalid_argument("The matrix has the wrong dimensions.");
    }
    this->matrix = std::move(matrix);
}

template <typename Derived>
const Derived &Operator<Derived>::derived() const {
    return static_cast<const Derived &>(*this);
}

template <typename Derived>
Derived &Operator<Derived>::derived_mutable() {
    return static_cast<Derived &>(*this);
}

template <typename Derived>
std::shared_ptr<const typename Operator<Derived>::basis_t> Operator<Derived>::get_basis() const {
    return basis;
}

template <typename Derived>
std::shared_ptr<const typename Operator<Derived>::basis_t> &Operator<Derived>::get_basis() {
    return basis;
}

template <typename Derived>
const Eigen::SparseMatrix<typename Operator<Derived>::scalar_t, Eigen::RowMajor> &
Operator<Derived>::get_matrix() const {
    return matrix;
}

template <typename Derived>
Eigen::SparseMatrix<typename Operator<Derived>::scalar_t, Eigen::RowMajor> &
Operator<Derived>::get_matrix() {
    return matrix;
}

template <typename Derived>
const Transformation<typename Operator<Derived>::scalar_t> &
Operator<Derived>::get_transformation() const {
    return basis->get_transformation();
}

template <typename Derived>
Transformation<typename Operator<Derived>::scalar_t>
Operator<Derived>::get_rotator(real_t alpha, real_t beta, real_t gamma) const {
    return basis->get_rotator(alpha, beta, gamma);
}

template <typename Derived>
Sorting Operator<Derived>::get_sorter(const std::vector<TransformationType> &labels) const {
    basis->perform_sorter_checks(labels);

    // Split labels into three parts (one before SORT_BY_ENERGY, one with SORT_BY_ENERGY, and one
    // after)
    auto it = std::find(labels.begin(), labels.end(), TransformationType::SORT_BY_ENERGY);
    std::vector<TransformationType> before_energy(labels.begin(), it);
    bool contains_energy = (it != labels.end());
    std::vector<TransformationType> after_energy(contains_energy ? it + 1 : labels.end(),
                                                 labels.end());

    // Initialize transformation
    Sorting transformation;
    transformation.matrix.resize(matrix.rows());
    transformation.matrix.setIdentity();

    // Apply sorting for labels before SORT_BY_ENERGY
    if (!before_energy.empty()) {
        basis->get_sorter_without_checks(before_energy, transformation);
    }

    // Apply SORT_BY_ENERGY if present
    if (contains_energy) {
        std::vector<real_t> energies_of_states;
        energies_of_states.reserve(matrix.rows());
        for (int i = 0; i < matrix.rows(); ++i) {
            energies_of_states.push_back(std::real(matrix.coeff(i, i)));
        }

        std::stable_sort(
            transformation.matrix.indices().data(),
            transformation.matrix.indices().data() + transformation.matrix.indices().size(),
            [&](int i, int j) { return energies_of_states[i] < energies_of_states[j]; });

        transformation.transformation_type.push_back(TransformationType::SORT_BY_ENERGY);
    }

    // Apply sorting for labels after SORT_BY_ENERGY
    if (!after_energy.empty()) {
        basis->get_sorter_without_checks(after_energy, transformation);
    }

    // Check if all labels have been used for sorting
    if (labels != transformation.transformation_type) {
        throw std::invalid_argument("The states could not be sorted by all the requested labels.");
    }

    return transformation;
}

template <typename Derived>
std::vector<IndicesOfBlock>
Operator<Derived>::get_indices_of_blocks(const std::vector<TransformationType> &labels) const {
    basis->perform_sorter_checks(labels);

    std::set<TransformationType> unique_labels(labels.begin(), labels.end());
    basis->perform_blocks_checks(unique_labels);

    // Split labels into two parts (one with SORT_BY_ENERGY and one without)
    auto it = unique_labels.find(TransformationType::SORT_BY_ENERGY);
    bool contains_energy = (it != unique_labels.end());
    if (contains_energy) {
        unique_labels.erase(it);
    }

    // Initialize blocks
    IndicesOfBlocksCreator blocks_creator({0, static_cast<size_t>(matrix.rows())});

    // Handle all labels except SORT_BY_ENERGY
    if (!unique_labels.empty()) {
        basis->get_indices_of_blocks_without_checks(unique_labels, blocks_creator);
    }

    // Handle SORT_BY_ENERGY if present
    if (contains_energy) {
        scalar_t last_energy = std::real(matrix.coeff(0, 0));
        for (int i = 0; i < matrix.rows(); ++i) {
            if (std::real(matrix.coeff(i, i)) != last_energy) {
                blocks_creator.add(i);
                last_energy = std::real(matrix.coeff(i, i));
            }
        }
    }

    return blocks_creator.create();
}

template <typename Derived>
Derived Operator<Derived>::transformed(
    const Transformation<typename Operator<Derived>::scalar_t> &transformation) const {
    auto transformed = derived();
    if (matrix.cols() == 0) {
        return transformed;
    }
    transformed.matrix = transformation.matrix.adjoint() * matrix * transformation.matrix;
    transformed.basis = basis->transformed(transformation);
    return transformed;
}

template <typename Derived>
Derived Operator<Derived>::transformed(const Sorting &transformation) const {
    auto transformed = derived();
    if (matrix.cols() == 0) {
        return transformed;
    }
    transformed.matrix = matrix.twistedBy(transformation.matrix.inverse());
    transformed.basis = basis->transformed(transformation);
    return transformed;
}

// Overloaded operators
template <typename Derived>
Derived operator*(const typename Operator<Derived>::scalar_t &lhs, const Operator<Derived> &rhs) {
    Derived result = rhs.derived();
    result.matrix *= lhs;
    return result;
}

template <typename Derived>
Derived operator*(const Operator<Derived> &lhs, const typename Operator<Derived>::scalar_t &rhs) {
    Derived result = lhs.derived();
    result.matrix *= rhs;
    return result;
}

template <typename Derived>
Derived operator/(const Operator<Derived> &lhs, const typename Operator<Derived>::scalar_t &rhs) {
    Derived result = lhs.derived();
    result.matrix /= rhs;
    return result;
}

template <typename Derived>
Derived &operator+=(Operator<Derived> &lhs, const Operator<Derived> &rhs) {
    if (lhs.basis != rhs.basis) {
        throw std::invalid_argument("The basis of the operators is not the same.");
    }
    lhs.matrix += rhs.matrix;
    return lhs.derived_mutable();
}

template <typename Derived>
Derived &operator-=(Operator<Derived> &lhs, const Operator<Derived> &rhs) {
    if (lhs.basis != rhs.basis) {
        throw std::invalid_argument("The basis of the operators is not the same.");
    }
    lhs.matrix -= rhs.matrix;
    return lhs.derived_mutable();
}

template <typename Derived>
Derived operator+(const Operator<Derived> &lhs, const Operator<Derived> &rhs) {
    if (lhs.basis != rhs.basis) {
        throw std::invalid_argument("The basis of the operators is not the same.");
    }
    Derived result = lhs.derived();
    result.matrix += rhs.matrix;
    return result;
}

template <typename Derived>
Derived operator-(const Operator<Derived> &lhs, const Operator<Derived> &rhs) {
    if (lhs.basis != rhs.basis) {
        throw std::invalid_argument("The basis of the operators is not the same.");
    }
    Derived result = lhs.derived();
    result.matrix -= rhs.matrix;
    return result;
}

// Explicit instantiations
// NOLINTBEGIN(bugprone-macro-parentheses, cppcoreguidelines-macro-usage)
#define INSTANTIATE_OPERATOR_HELPER(SCALAR, TYPE)                                                  \
    template class Operator<TYPE<SCALAR>>;                                                         \
    template TYPE<SCALAR> operator*(const SCALAR &lhs, const Operator<TYPE<SCALAR>> &rhs);         \
    template TYPE<SCALAR> operator*(const Operator<TYPE<SCALAR>> &lhs, const SCALAR &rhs);         \
    template TYPE<SCALAR> operator/(const Operator<TYPE<SCALAR>> &lhs, const SCALAR &rhs);         \
    template TYPE<SCALAR> &operator+=(Operator<TYPE<SCALAR>> &lhs,                                 \
                                      const Operator<TYPE<SCALAR>> &rhs);                          \
    template TYPE<SCALAR> &operator-=(Operator<TYPE<SCALAR>> &lhs,                                 \
                                      const Operator<TYPE<SCALAR>> &rhs);                          \
    template TYPE<SCALAR> operator+(const Operator<TYPE<SCALAR>> &lhs,                             \
                                    const Operator<TYPE<SCALAR>> &rhs);                            \
    template TYPE<SCALAR> operator-(const Operator<TYPE<SCALAR>> &lhs,                             \
                                    const Operator<TYPE<SCALAR>> &rhs);
#define INSTANTIATE_OPERATOR(SCALAR)                                                               \
    INSTANTIATE_OPERATOR_HELPER(SCALAR, OperatorAtom)                                              \
    INSTANTIATE_OPERATOR_HELPER(SCALAR, OperatorPair)
// NOLINTEND(bugprone-macro-parentheses, cppcoreguidelines-macro-usage)

INSTANTIATE_OPERATOR(double)
INSTANTIATE_OPERATOR(std::complex<double>)

#undef INSTANTIATE_OPERATOR_HELPER
#undef INSTANTIATE_OPERATOR

} // namespace pairinteraction
