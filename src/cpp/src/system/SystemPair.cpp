// SPDX-FileCopyrightText: 2024 Pairinteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "pairinteraction/system/SystemPair.hpp"

#include "pairinteraction/basis/BasisAtom.hpp"
#include "pairinteraction/basis/BasisPair.hpp"
#include "pairinteraction/enums/OperatorType.hpp"
#include "pairinteraction/enums/Parity.hpp"
#include "pairinteraction/enums/TransformationType.hpp"
#include "pairinteraction/ket/KetAtom.hpp"
#include "pairinteraction/ket/KetPair.hpp"
#include "pairinteraction/operator/OperatorAtom.hpp"
#include "pairinteraction/operator/OperatorPair.hpp"
#include "pairinteraction/system/GreenTensor.hpp"
#include "pairinteraction/system/SystemAtom.hpp"
#include "pairinteraction/utils/Range.hpp"
#include "pairinteraction/utils/eigen_assertion.hpp"
#include "pairinteraction/utils/eigen_compat.hpp"
#include "pairinteraction/utils/spherical.hpp"
#include "pairinteraction/utils/streamed.hpp"
#include "pairinteraction/utils/tensor.hpp"
#include "pairinteraction/utils/traits.hpp"

#include <Eigen/SparseCore>
#include <array>
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
};

template <typename Scalar>
GreenTensor<Scalar> construct_green_tensor(
    const std::array<typename traits::NumTraits<Scalar>::real_t, 3> &distance_vector, int order) {
    // https://doi.org/10.1103/PhysRevA.96.062509
    // https://doi.org/10.1103/PhysRevA.82.010901
    // https://en.wikipedia.org/wiki/Table_of_spherical_harmonics

    using real_t = typename traits::NumTraits<Scalar>::real_t;

    GreenTensor<Scalar> green_tensor;

    // Normalize the distance vector, return zero green tensor if the distance is infinity
    Eigen::Map<const Eigen::Vector3<real_t>> vector_map(distance_vector.data(),
                                                        distance_vector.size());
    real_t distance = vector_map.norm();
    SPDLOG_DEBUG("Interatomic distance: {}", distance);
    if (distance == std::numeric_limits<real_t>::infinity()) {
        return green_tensor;
    }
    Eigen::Vector3<real_t> vector_normalized = vector_map / distance;

    // Dyadic green function of dipole-dipole interaction
    if (order >= 3) {
        Eigen::Matrix3<Scalar> entries = Eigen::Matrix3<real_t>::Identity() -
            3 * vector_normalized * vector_normalized.transpose();
        entries /= std::pow(distance, 3);

        green_tensor.set_entries(1, 1, entries.template cast<Scalar>());
    }

    // Dyadic green function of dipole-quadrupole interaction
    if (order >= 4) {
        Eigen::Matrix<real_t, 3, 9> entries = Eigen::Matrix<real_t, 3, 9>::Zero();
        for (Eigen::Index q = 0; q < 3; ++q) {
            Eigen::Index row = q;
            for (Eigen::Index j = 0; j < 3; ++j) {
                for (Eigen::Index i = 0; i < 3; ++i) {
                    Eigen::Index col = 3 * j + i;
                    entries(row, col) +=
                        15 * vector_normalized[q] * vector_normalized[j] * vector_normalized[i];
                    if (i == j) {
                        entries(row, col) += -3 * vector_normalized[q];
                    }
                    if (i == q) {
                        entries(row, col) += -3 * vector_normalized[j];
                    }
                    if (j == q) {
                        entries(row, col) += -3 * vector_normalized[i];
                    }
                }
            }
        }
        entries /= std::pow(distance, 4);

        green_tensor.set_entries(1, 2, entries.template cast<Scalar>());
    }

    // Dyadic green function of quadrupole-dipole interaction
    if (order >= 4) {
        Eigen::Matrix<real_t, 9, 3> entries = Eigen::Matrix<real_t, 9, 3>::Zero();
        for (Eigen::Index q = 0; q < 3; ++q) {
            for (Eigen::Index j = 0; j < 3; ++j) {
                Eigen::Index row = 3 * q + j;
                for (Eigen::Index i = 0; i < 3; ++i) {
                    Eigen::Index col = i;
                    entries(row, col) +=
                        -15 * vector_normalized[q] * vector_normalized[j] * vector_normalized[i];
                    if (i == j) {
                        entries(row, col) += 3 * vector_normalized[q];
                    }
                    if (i == q) {
                        entries(row, col) += 3 * vector_normalized[j];
                    }
                    if (j == q) {
                        entries(row, col) += 3 * vector_normalized[i];
                    }
                }
            }
        }
        entries /= std::pow(distance, 4);

        green_tensor.set_entries(2, 1, entries.template cast<Scalar>());
    }

    // Dyadic green function of quadrupole-quadrupole interaction
    if (order >= 5) {
        SPDLOG_WARN("Quadrupole-quadrupole interaction is considered but "
                    "not dipole-octupole interaction although this interaction would be "
                    "of the same order. We plan to implement dipole-octupole interaction "
                    "in the future.");

        Eigen::Matrix<real_t, 9, 9> entries = Eigen::Matrix<real_t, 9, 9>::Zero();
        for (Eigen::Index q = 0; q < 3; ++q) {
            for (Eigen::Index j = 0; j < 3; ++j) {
                Eigen::Index row = 3 * q + j;
                for (Eigen::Index i = 0; i < 3; ++i) {
                    for (Eigen::Index k = 0; k < 3; ++k) {
                        Eigen::Index col = 3 * i + k;
                        entries(row, col) += 105 * vector_normalized[q] * vector_normalized[j] *
                            vector_normalized[i] * vector_normalized[k];
                        if (i == j) {
                            entries(row, col) += -15 * vector_normalized[q] * vector_normalized[k];
                        }
                        if (i == q) {
                            entries(row, col) += -15 * vector_normalized[j] * vector_normalized[k];
                        }
                        if (j == q) {
                            entries(row, col) += -15 * vector_normalized[i] * vector_normalized[k];
                        }
                        if (k == q) {
                            entries(row, col) += -15 * vector_normalized[j] * vector_normalized[i];
                        }
                        if (k == j) {
                            entries(row, col) += -15 * vector_normalized[q] * vector_normalized[i];
                        }
                        if (k == i) {
                            entries(row, col) += -15 * vector_normalized[q] * vector_normalized[j];
                        }
                        if (q == k && i == j) {
                            entries(row, col) += 3;
                        }
                        if (i == k && j == q) {
                            entries(row, col) += 3;
                        }
                        if (j == k && i == q) {
                            entries(row, col) += 3;
                        }
                    }
                }
            }
        }
        entries /= std::pow(distance, 5);

        green_tensor.set_entries(2, 2, entries.template cast<Scalar>());
    }

    return green_tensor;
}

template <typename Scalar>
OperatorMatrices<Scalar>
construct_operator_matrices(const GreenTensor<Scalar> &green_tensor,
                            const std::shared_ptr<const BasisAtom<Scalar>> &basis1,
                            const std::shared_ptr<const BasisAtom<Scalar>> &basis2) {
    OperatorMatrices<Scalar> op;

    if (!green_tensor.get_entries(1, 1).empty() || !green_tensor.get_entries(1, 2).empty()) {
        op.d1.push_back(
            -OperatorAtom<Scalar>(basis1, OperatorType::ELECTRIC_DIPOLE, 1).get_matrix());
        op.d1.push_back(
            OperatorAtom<Scalar>(basis1, OperatorType::ELECTRIC_DIPOLE, 0).get_matrix());
        op.d1.push_back(
            -OperatorAtom<Scalar>(basis1, OperatorType::ELECTRIC_DIPOLE, -1).get_matrix());
    }

    if (!green_tensor.get_entries(1, 1).empty() || !green_tensor.get_entries(2, 1).empty()) {
        op.d2.push_back(
            OperatorAtom<Scalar>(basis2, OperatorType::ELECTRIC_DIPOLE, -1).get_matrix());
        op.d2.push_back(
            OperatorAtom<Scalar>(basis2, OperatorType::ELECTRIC_DIPOLE, 0).get_matrix());
        op.d2.push_back(
            OperatorAtom<Scalar>(basis2, OperatorType::ELECTRIC_DIPOLE, 1).get_matrix());
    }

    if (!green_tensor.get_entries(2, 2).empty() || !green_tensor.get_entries(2, 1).empty()) {
        op.q1.push_back(
            OperatorAtom<Scalar>(basis1, OperatorType::ELECTRIC_QUADRUPOLE, 2).get_matrix());
        op.q1.push_back(
            -OperatorAtom<Scalar>(basis1, OperatorType::ELECTRIC_QUADRUPOLE, 1).get_matrix());
        op.q1.push_back(
            OperatorAtom<Scalar>(basis1, OperatorType::ELECTRIC_QUADRUPOLE, 0).get_matrix());
        op.q1.push_back(
            -OperatorAtom<Scalar>(basis1, OperatorType::ELECTRIC_QUADRUPOLE, -1).get_matrix());
        op.q1.push_back(
            OperatorAtom<Scalar>(basis1, OperatorType::ELECTRIC_QUADRUPOLE, -2).get_matrix());
        op.q1.push_back(
            OperatorAtom<Scalar>(basis1, OperatorType::ELECTRIC_QUADRUPOLE_ZERO, 0).get_matrix());
    }

    if (!green_tensor.get_entries(2, 2).empty() || !green_tensor.get_entries(1, 2).empty()) {
        op.q2.push_back(
            OperatorAtom<Scalar>(basis2, OperatorType::ELECTRIC_QUADRUPOLE, -2).get_matrix());
        op.q2.push_back(
            OperatorAtom<Scalar>(basis2, OperatorType::ELECTRIC_QUADRUPOLE, -1).get_matrix());
        op.q2.push_back(
            OperatorAtom<Scalar>(basis2, OperatorType::ELECTRIC_QUADRUPOLE, 0).get_matrix());
        op.q2.push_back(
            OperatorAtom<Scalar>(basis2, OperatorType::ELECTRIC_QUADRUPOLE, 1).get_matrix());
        op.q2.push_back(
            OperatorAtom<Scalar>(basis2, OperatorType::ELECTRIC_QUADRUPOLE, 2).get_matrix());
        op.q2.push_back(
            OperatorAtom<Scalar>(basis2, OperatorType::ELECTRIC_QUADRUPOLE_ZERO, 0).get_matrix());
    }

    return op;
}

template <typename Scalar>
SystemPair<Scalar>::SystemPair(std::shared_ptr<const basis_t> basis)
    : System<SystemPair<Scalar>>(std::move(basis)) {}

template <typename Scalar>
SystemPair<Scalar> &SystemPair<Scalar>::set_order(int value) {
    this->hamiltonian_requires_construction = true;
    if (value < 3 || value > 5) {
        throw std::invalid_argument("The order must be 3, 4, or 5.");
    }
    order = value;
    return *this;
}

template <typename Scalar>
SystemPair<Scalar> &SystemPair<Scalar>::set_distance_vector(const std::array<real_t, 3> &vector) {
    this->hamiltonian_requires_construction = true;

    constexpr real_t numerical_precision = 100 * std::numeric_limits<real_t>::epsilon();
    if (!traits::NumTraits<Scalar>::is_complex_v &&
        std::abs(distance_vector[1]) > numerical_precision) {
        throw std::invalid_argument(
            "The distance vector must not have a y-component if the scalar type is real.");
    }

    distance_vector = vector;

    return *this;
}

template <typename Scalar>
void SystemPair<Scalar>::construct_hamiltonian() const {
    auto basis = this->hamiltonian->get_basis();
    auto basis1 = basis->get_basis1();
    auto basis2 = basis->get_basis2();

    auto green_tensor = construct_green_tensor<Scalar>(distance_vector, order);
    auto op = construct_operator_matrices(green_tensor, basis1, basis2);

    // Construct the unperturbed Hamiltonian
    this->hamiltonian = std::make_unique<OperatorPair<Scalar>>(basis, OperatorType::ENERGY);
    this->hamiltonian_is_diagonal = true;
    bool sort_by_quantum_number_f = basis->has_quantum_number_f();
    bool sort_by_quantum_number_m = basis->has_quantum_number_m();
    bool sort_by_parity = basis->has_parity();

    // Dipole-dipole interaction
    for (const auto &entry : green_tensor.get_entries(1, 1)) {
        if (auto ce = std::get_if<typename GreenTensor<Scalar>::ConstantEntry>(&entry)) {
            this->hamiltonian->get_matrix() += ce->val() *
                utils::calculate_tensor_product(basis, basis, op.d1[ce->row()], op.d2[ce->col()]);
            if (ce->row() != ce->col()) {
                sort_by_quantum_number_m = false;
            }
            this->hamiltonian_is_diagonal = false;
            sort_by_quantum_number_f = false;
        } else {
            throw std::invalid_argument("Omega-dependent entries are not yet supported.");
        }
    }

    // Dipole-quadrupole interaction
    for (const auto &entry : green_tensor.get_entries(1, 2)) {
        if (auto ce = std::get_if<typename GreenTensor<Scalar>::ConstantEntry>(&entry)) {
            this->hamiltonian->get_matrix() += ce->val() *
                utils::calculate_tensor_product(basis, basis, op.d1[ce->row()], op.q2[ce->col()]);
            if (ce->row() != ce->col() - 1) {
                sort_by_quantum_number_m = false;
            }
            this->hamiltonian_is_diagonal = false;
            sort_by_quantum_number_f = false;
        } else {
            throw std::invalid_argument("Omega-dependent entries are not yet supported.");
        }
    }

    // Quadrupole-dipole interaction
    for (const auto &entry : green_tensor.get_entries(2, 1)) {
        if (auto ce = std::get_if<typename GreenTensor<Scalar>::ConstantEntry>(&entry)) {
            this->hamiltonian->get_matrix() += ce->val() *
                utils::calculate_tensor_product(basis, basis, op.q1[ce->row()], op.d2[ce->col()]);
            if (ce->row() - 1 != ce->col()) {
                sort_by_quantum_number_m = false;
            }
            this->hamiltonian_is_diagonal = false;
            sort_by_quantum_number_f = false;
        } else {
            throw std::invalid_argument("Omega-dependent entries are not yet supported.");
        }
    }

    // Quadrupole-quadrupole interaction
    for (const auto &entry : green_tensor.get_entries(2, 2)) {
        if (auto ce = std::get_if<typename GreenTensor<Scalar>::ConstantEntry>(&entry)) {
            this->hamiltonian->get_matrix() += ce->val() *
                utils::calculate_tensor_product(basis, basis, op.q1[ce->row()], op.q2[ce->col()]);
            if (ce->row() != ce->col()) {
                sort_by_quantum_number_m = false;
            }
            this->hamiltonian_is_diagonal = false;
            sort_by_quantum_number_f = false;
        } else {
            throw std::invalid_argument("Omega-dependent entries are not yet supported.");
        }
    }

    // Store which labels can be used to block-diagonalize the Hamiltonian
    this->blockdiagonalizing_labels.clear();
    if (sort_by_quantum_number_f) {
        this->blockdiagonalizing_labels.push_back(TransformationType::SORT_BY_QUANTUM_NUMBER_F);
    }
    if (sort_by_quantum_number_m) {
        this->blockdiagonalizing_labels.push_back(TransformationType::SORT_BY_QUANTUM_NUMBER_M);
    }
    if (sort_by_parity) {
        this->blockdiagonalizing_labels.push_back(TransformationType::SORT_BY_PARITY);
    }
}

// Explicit instantiations
template class SystemPair<double>;
template class SystemPair<std::complex<double>>;
} // namespace pairinteraction
