// SPDX-FileCopyrightText: 2024 PairInteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "pairinteraction/system/SystemPair.hpp"

#include "pairinteraction/basis/BasisAtom.hpp"
#include "pairinteraction/basis/BasisPair.hpp"
#include "pairinteraction/database/Database.hpp"
#include "pairinteraction/enums/OperatorType.hpp"
#include "pairinteraction/enums/Parity.hpp"
#include "pairinteraction/enums/SorterType.hpp"
#include "pairinteraction/ket/KetAtom.hpp"
#include "pairinteraction/ket/KetPair.hpp"
#include "pairinteraction/system/GreenTensorInterpolator.hpp"
#include "pairinteraction/system/SystemAtom.hpp"
#include "pairinteraction/utils/Range.hpp"
#include "pairinteraction/utils/eigen_assertion.hpp"
#include "pairinteraction/utils/eigen_compat.hpp"
#include "pairinteraction/utils/operator.hpp"
#include "pairinteraction/utils/streamed.hpp"
#include "pairinteraction/utils/tensor.hpp"
#include "pairinteraction/utils/traits.hpp"

#include <Eigen/SparseCore>
#include <array>
#include <cmath>
#include <complex>
#include <limits>
#include <memory>
#include <spdlog/spdlog.h>
#include <vector>

namespace pairinteraction {
template <typename Scalar>
struct OperatorMatrices {
    std::vector<Eigen::SparseMatrix<Scalar, Eigen::RowMajor>> d1;
    std::vector<Eigen::SparseMatrix<Scalar, Eigen::RowMajor>> d2;
    std::vector<Eigen::SparseMatrix<Scalar, Eigen::RowMajor>> q1;
    std::vector<Eigen::SparseMatrix<Scalar, Eigen::RowMajor>> q2;
    std::vector<Eigen::SparseMatrix<Scalar, Eigen::RowMajor>> o1;
    std::vector<Eigen::SparseMatrix<Scalar, Eigen::RowMajor>> o2;
};

template <typename Scalar>
OperatorMatrices<Scalar>
construct_operator_matrices(const GreenTensorInterpolator<Scalar> &green_tensor_interpolator,
                            const std::shared_ptr<const BasisAtom<Scalar>> &basis1,
                            const std::shared_ptr<const BasisAtom<Scalar>> &basis2) {
    // Helper function for constructing matrices of spherical harmonics operators
    auto get_matrices = [](auto basis, OperatorType type, std::initializer_list<int> m,
                           bool conjugate) {
        std::vector<Eigen::SparseMatrix<Scalar, Eigen::RowMajor>> matrices;
        matrices.reserve(m.size());
        int factor = conjugate ? -1 : 1;
        std::transform(m.begin(), m.end(), std::back_inserter(matrices), [&](int q) {
            auto matrix_elements = (std::pow(factor, q) *
                                    basis->get_database().get_matrix_elements_in_canonical_basis(
                                        basis, basis, type, factor * q))
                                       .eval();
            return (basis->get_coefficients().adjoint() * matrix_elements *
                    basis->get_coefficients())
                .eval();
        });
        return matrices;
    };

    OperatorMatrices<Scalar> op;

    // Operator matrices for Rydberg-Rydberg interaction
    if (!green_tensor_interpolator.get_spherical_entries(1, 1).empty() ||
        !green_tensor_interpolator.get_spherical_entries(1, 2).empty() ||
        !green_tensor_interpolator.get_spherical_entries(1, 3).empty()) {
        op.d1 = get_matrices(basis1, OperatorType::ELECTRIC_DIPOLE, {-1, 0, +1}, true);
    }
    if (!green_tensor_interpolator.get_spherical_entries(1, 1).empty() ||
        !green_tensor_interpolator.get_spherical_entries(2, 1).empty() ||
        !green_tensor_interpolator.get_spherical_entries(3, 1).empty()) {
        op.d2 = get_matrices(basis2, OperatorType::ELECTRIC_DIPOLE, {-1, 0, +1}, false);
    }
    if (!green_tensor_interpolator.get_spherical_entries(2, 2).empty() ||
        !green_tensor_interpolator.get_spherical_entries(2, 1).empty()) {
        op.q1 = get_matrices(basis1, OperatorType::ELECTRIC_QUADRUPOLE, {-2, -1, 0, +1, +2}, true);
        op.q1.push_back(get_matrices(basis1, OperatorType::ELECTRIC_QUADRUPOLE_ZERO, {0}, true)[0]);
    }
    if (!green_tensor_interpolator.get_spherical_entries(2, 2).empty() ||
        !green_tensor_interpolator.get_spherical_entries(1, 2).empty()) {
        op.q2 = get_matrices(basis2, OperatorType::ELECTRIC_QUADRUPOLE, {-2, -1, 0, +1, +2}, false);
        op.q2.push_back(
            get_matrices(basis2, OperatorType::ELECTRIC_QUADRUPOLE_ZERO, {0}, false)[0]);
    }
    // In contrast to the quadrupole operators, no trace operator must be appended because the
    // cartesian-to-spherical transformator for kappa == 3 does not contain trace rows.
    if (!green_tensor_interpolator.get_spherical_entries(3, 1).empty()) {
        op.o1 = get_matrices(basis1, OperatorType::ELECTRIC_OCTUPOLE, {-3, -2, -1, 0, +1, +2, +3},
                             true);
    }
    if (!green_tensor_interpolator.get_spherical_entries(1, 3).empty()) {
        op.o2 = get_matrices(basis2, OperatorType::ELECTRIC_OCTUPOLE, {-3, -2, -1, 0, +1, +2, +3},
                             false);
    }

    return op;
}

template <typename Scalar>
SystemPair<Scalar>::SystemPair(std::shared_ptr<const basis_t> basis)
    : System<SystemPair<Scalar>>(std::move(basis)) {}

template <typename Scalar>
SystemPair<Scalar> &SystemPair<Scalar>::set_interaction_order(int value) {
    this->hamiltonian_requires_construction = true;

    if (value < 3 || value > 5) {
        throw std::invalid_argument("The order must be 3, 4, or 5.");
    }

    if (green_tensor_interpolator) {
        throw std::invalid_argument(
            "Cannot set interaction order if a user-defined green tensor interpolator is set.");
    }

    interaction_order = value;

    return *this;
}

template <typename Scalar>
SystemPair<Scalar> &SystemPair<Scalar>::set_distance_vector(const std::array<real_t, 3> &vector) {
    this->hamiltonian_requires_construction = true;

    if (!traits::NumTraits<Scalar>::is_complex_v && vector[1] != 0) {
        throw std::invalid_argument(
            "The distance vector must not have a y-component if the scalar type is real.");
    }

    if (green_tensor_interpolator) {
        throw std::invalid_argument(
            "Cannot set distance vector if a user-defined green tensor interpolator is set.");
    }

    distance_vector = vector;

    return *this;
}

template <typename Scalar>
SystemPair<Scalar> &SystemPair<Scalar>::set_green_tensor_interpolator(
    const std::shared_ptr<const GreenTensorInterpolator<Scalar>> &green_tensor_interpolator) {
    this->hamiltonian_requires_construction = true;

    if (std::isfinite(distance_vector[0]) && std::isfinite(distance_vector[1]) &&
        std::isfinite(distance_vector[2])) {
        throw std::invalid_argument(
            "Cannot set green tensor interpolator if a finite distance vector is set.");
    }

    this->green_tensor_interpolator = green_tensor_interpolator;

    return *this;
}

template <typename Scalar>
void SystemPair<Scalar>::construct_hamiltonian() const {
    auto basis1 = this->basis->get_basis1();
    auto basis2 = this->basis->get_basis2();

    std::shared_ptr<const GreenTensorInterpolator<Scalar>> green_tensor_interpolator_ptr;
    if (green_tensor_interpolator) {
        green_tensor_interpolator_ptr = green_tensor_interpolator;
    } else {
        green_tensor_interpolator_ptr = std::make_shared<const GreenTensorInterpolator<Scalar>>(
            GreenTensorInterpolator<Scalar>::from_multipole_expansion(distance_vector,
                                                                      interaction_order));
    }

    auto op = construct_operator_matrices(*green_tensor_interpolator_ptr, basis1, basis2);

    // Construct the unperturbed Hamiltonian in the canonical pair basis
    this->matrix = utils::get_energies_in_canonical_basis(this->basis);

    this->hamiltonian_is_diagonal = false;
    bool sort_by_quantum_number_f = this->basis->has_quantum_number_f();
    bool sort_by_quantum_number_m = this->basis->has_quantum_number_m();
    bool sort_by_parity = this->basis->has_parity();

    // Add Rydberg-Rydberg interaction via Green tensor
    // H_RR = Σ_{ij} D_1,left[i] * G_{ij} * D_2,right[j]
    // where D_1,left uses conjugated convention and
    // D_2,right uses normal convention.

    // Helper function for adding Rydberg-Rydberg interaction.
    auto add_interaction = [this, &green_tensor_interpolator_ptr, &sort_by_quantum_number_f,
                            &sort_by_quantum_number_m, &sort_by_parity](
                               const auto &op1, const auto &op2, int kappa1, int kappa2) {
        const auto &entries = green_tensor_interpolator_ptr->get_spherical_entries(kappa1, kappa2);

        for (const auto &entry : entries) {
            if (std::holds_alternative<
                    typename GreenTensorInterpolator<Scalar>::OmegaDependentEntry>(entry)) {
                throw std::logic_error(
                    "Green tensor with omega dependent entries is currently not supported.");
            }

            const auto &constant_entry =
                std::get<typename GreenTensorInterpolator<Scalar>::ConstantEntry>(entry);
            this->matrix += constant_entry.val() *
                utils::calculate_tensor_product_in_canonical_basis(this->basis, this->basis,
                                                                   op1[constant_entry.row()],
                                                                   op2[constant_entry.col()]);

            sort_by_quantum_number_f = false;
            if (constant_entry.row() != constant_entry.col() + kappa1 - kappa2) {
                sort_by_quantum_number_m = false;
            }
            if ((kappa1 + kappa2) % 2 != 0) {
                sort_by_parity = false;
            }
        }
    };

    // Dipole-dipole interaction
    add_interaction(op.d1, op.d2, 1, 1);

    // Dipole-quadrupole interaction
    add_interaction(op.d1, op.q2, 1, 2);

    // Quadrupole-dipole interaction
    add_interaction(op.q1, op.d2, 2, 1);

    // Quadrupole-quadrupole interaction
    add_interaction(op.q1, op.q2, 2, 2);

    // Dipole-octupole interaction
    add_interaction(op.d1, op.o2, 1, 3);

    // Octupole-dipole interaction
    add_interaction(op.o1, op.d2, 3, 1);

    // Transform from the canonical basis into the actual basis
    this->matrix =
        this->basis->get_coefficients().adjoint() * this->matrix * this->basis->get_coefficients();

    // Store which labels can be used to block-diagonalize the Hamiltonian
    this->blockdiagonalizing_labels.clear();
    if (sort_by_quantum_number_f) {
        this->blockdiagonalizing_labels.push_back(SorterType::QUANTUM_NUMBER_F);
    }
    if (sort_by_quantum_number_m) {
        this->blockdiagonalizing_labels.push_back(SorterType::QUANTUM_NUMBER_M);
    }
    if (sort_by_parity) {
        this->blockdiagonalizing_labels.push_back(SorterType::PARITY);
    }
}

// Explicit instantiations
template class SystemPair<double>;
template class SystemPair<std::complex<double>>;
} // namespace pairinteraction
