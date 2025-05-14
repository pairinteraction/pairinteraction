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
    const std::array<typename traits::NumTraits<Scalar>::real_t, 3> &distance_vector,
    int interaction_order) {
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
    if (!std::isfinite(distance)) {
        return green_tensor;
    }
    Eigen::Vector3<real_t> unitvec = vector_map / distance;

    // Dyadic green function of dipole-dipole interaction
    if (interaction_order >= 3) {
        Eigen::Matrix3<Scalar> entries =
            Eigen::Matrix3<real_t>::Identity() - 3 * unitvec * unitvec.transpose();

        green_tensor.set_entries(1, 1, (entries / std::pow(distance, 3)).template cast<Scalar>());
    }

    // Dyadic green function of dipole-quadrupole interaction
    if (interaction_order >= 4) {
        Eigen::Matrix<real_t, 3, 9> entries = Eigen::Matrix<real_t, 3, 9>::Zero();
        for (Eigen::Index q = 0; q < 3; ++q) {
            Eigen::Index row = q;
            for (Eigen::Index j = 0; j < 3; ++j) {
                for (Eigen::Index i = 0; i < 3; ++i) {
                    Eigen::Index col = 3 * j + i;
                    real_t v = 15 * unitvec[q] * unitvec[j] * unitvec[i];
                    if (i == j) v += -3 * unitvec[q];
                    if (i == q) v += -3 * unitvec[j];
                    if (j == q) v += -3 * unitvec[i];
                    entries(row, col) += v;
                }
            }
        }

        green_tensor.set_entries(1, 2, (entries / std::pow(distance, 4)).template cast<Scalar>());
    }

    // Dyadic green function of quadrupole-dipole interaction
    if (interaction_order >= 4) {
        Eigen::Matrix<real_t, 9, 3> entries = Eigen::Matrix<real_t, 9, 3>::Zero();
        for (Eigen::Index q = 0; q < 3; ++q) {
            for (Eigen::Index j = 0; j < 3; ++j) {
                Eigen::Index row = 3 * q + j;
                for (Eigen::Index i = 0; i < 3; ++i) {
                    Eigen::Index col = i;
                    real_t v = -15 * unitvec[q] * unitvec[j] * unitvec[i];
                    if (i == j) v += 3 * unitvec[q];
                    if (i == q) v += 3 * unitvec[j];
                    if (j == q) v += 3 * unitvec[i];
                    entries(row, col) += v;
                }
            }
        }

        green_tensor.set_entries(2, 1, (entries / std::pow(distance, 4)).template cast<Scalar>());
    }

    // Dyadic green function of quadrupole-quadrupole interaction
    if (interaction_order >= 5) {
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
                        real_t v = 105 * unitvec[q] * unitvec[j] * unitvec[i] * unitvec[k];
                        if (i == j) v += -15 * unitvec[q] * unitvec[k];
                        if (i == q) v += -15 * unitvec[j] * unitvec[k];
                        if (j == q) v += -15 * unitvec[i] * unitvec[k];
                        if (k == q) v += -15 * unitvec[j] * unitvec[i];
                        if (k == j) v += -15 * unitvec[q] * unitvec[i];
                        if (k == i) v += -15 * unitvec[q] * unitvec[j];
                        if (q == k && i == j) v += 3;
                        if (i == k && j == q) v += 3;
                        if (j == k && i == q) v += 3;
                        entries(row, col) += v;
                    }
                }
            }
        }

        green_tensor.set_entries(2, 2, (entries / std::pow(distance, 5)).template cast<Scalar>());
    }

    return green_tensor;
}

template <typename Scalar>
OperatorMatrices<Scalar>
construct_operator_matrices(const GreenTensor<Scalar> &green_tensor,
                            const std::shared_ptr<const BasisAtom<Scalar>> &basis1,
                            const std::shared_ptr<const BasisAtom<Scalar>> &basis2) {
    // Helper function for constructing matrices of spherical harmonics operators
    auto get_matrices = [](auto basis, OperatorType type, std::initializer_list<int> m,
                           bool conjugate) {
        std::vector<Eigen::SparseMatrix<Scalar, Eigen::RowMajor>> matrices;
        matrices.reserve(m.size());
        int factor = conjugate ? -1 : 1;
        std::transform(m.begin(), m.end(), std::back_inserter(matrices), [&](int q) {
            return (std::pow(factor, q) *
                    OperatorAtom<Scalar>(basis, type, factor * q).get_matrix())
                .eval();
        });
        return matrices;
    };

    OperatorMatrices<Scalar> op;

    // Operator matrices for Rydberg-Rydberg interaction
    if (!green_tensor.get_entries(1, 1).empty() || !green_tensor.get_entries(1, 2).empty()) {
        op.d1 = get_matrices(basis1, OperatorType::ELECTRIC_DIPOLE, {-1, 0, +1}, true);
    }
    if (!green_tensor.get_entries(1, 1).empty() || !green_tensor.get_entries(2, 1).empty()) {
        op.d2 = get_matrices(basis2, OperatorType::ELECTRIC_DIPOLE, {-1, 0, +1}, false);
    }
    if (!green_tensor.get_entries(2, 2).empty() || !green_tensor.get_entries(2, 1).empty()) {
        op.q1 = get_matrices(basis1, OperatorType::ELECTRIC_QUADRUPOLE, {-2, -1, 0, +1, +2}, true);
        op.q1.push_back(get_matrices(basis1, OperatorType::ELECTRIC_QUADRUPOLE_ZERO, {0}, true)[0]);
    }
    if (!green_tensor.get_entries(2, 2).empty() || !green_tensor.get_entries(1, 2).empty()) {
        op.q2 = get_matrices(basis2, OperatorType::ELECTRIC_QUADRUPOLE, {-2, -1, 0, +1, +2}, false);
        op.q2.push_back(
            get_matrices(basis2, OperatorType::ELECTRIC_QUADRUPOLE_ZERO, {0}, false)[0]);
    }

    return op;
}

// "overloaded" pattern for std::visit
template <class... Ts>
struct overloaded : Ts... {
    using Ts::operator()...;
};
template <class... Ts>
overloaded(Ts...) -> overloaded<Ts...>;

template <typename Scalar>
SystemPair<Scalar>::SystemPair(std::shared_ptr<const basis_t> basis)
    : System<SystemPair<Scalar>>(std::move(basis)) {}

template <typename Scalar>
SystemPair<Scalar> &SystemPair<Scalar>::set_interaction_order(int value) {
    this->hamiltonian_requires_construction = true;

    if (value < 3 || value > 5) {
        throw std::invalid_argument("The order must be 3, 4, or 5.");
    }

    if (user_defined_green_tensor) {
        throw std::invalid_argument(
            "Cannot set interaction order if a user-defined green tensor is set.");
    }

    interaction_order = value;

    return *this;
}

template <typename Scalar>
SystemPair<Scalar> &SystemPair<Scalar>::set_distance_vector(const std::array<real_t, 3> &vector) {
    this->hamiltonian_requires_construction = true;

    if (!traits::NumTraits<Scalar>::is_complex_v && distance_vector[1] != 0) {
        throw std::invalid_argument(
            "The distance vector must not have a y-component if the scalar type is real.");
    }

    if (user_defined_green_tensor) {
        throw std::invalid_argument(
            "Cannot set distance vector if a user-defined green tensor is set.");
    }

    distance_vector = vector;

    return *this;
}

template <typename Scalar>
SystemPair<Scalar> &
SystemPair<Scalar>::set_green_tensor(std::shared_ptr<const GreenTensor<Scalar>> &green_tensor) {
    this->hamiltonian_requires_construction = true;

    if (std::isfinite(distance_vector[0]) && std::isfinite(distance_vector[1]) &&
        std::isfinite(distance_vector[2])) {
        throw std::invalid_argument("Cannot set green tensor if a finite distance vector is set.");
    }

    user_defined_green_tensor = green_tensor;

    return *this;
}

template <typename Scalar>
void SystemPair<Scalar>::construct_hamiltonian() const {
    auto basis = this->hamiltonian->get_basis();
    auto basis1 = basis->get_basis1();
    auto basis2 = basis->get_basis2();

    std::shared_ptr<const GreenTensor<Scalar>> green_tensor_ptr;
    if (user_defined_green_tensor) {
        green_tensor_ptr = user_defined_green_tensor;
    } else {
        green_tensor_ptr = std::make_shared<const GreenTensor<Scalar>>(
            construct_green_tensor<Scalar>(distance_vector, interaction_order));
    }

    auto op = construct_operator_matrices(*green_tensor_ptr, basis1, basis2);

    // Construct the unperturbed Hamiltonian
    this->hamiltonian = std::make_unique<OperatorPair<Scalar>>(basis, OperatorType::ENERGY);
    this->hamiltonian_is_diagonal = true;
    bool sort_by_quantum_number_f = basis->has_quantum_number_f();
    bool sort_by_quantum_number_m = basis->has_quantum_number_m();
    bool sort_by_parity = basis->has_parity();

    // Store the energies (they are needed in case of Rydberg-Rydberg interaction with an
    // OmegaDependentEntry)
    auto energies = this->hamiltonian->get_matrix().diagonal().real();

    // Helper function for adding Rydberg-Rydberg interaction
    auto add_interaction = [this, &basis, &energies, &sort_by_quantum_number_f,
                            &sort_by_quantum_number_m](const auto &entries, const auto &op1,
                                                       const auto &op2, int delta) {
        for (const auto &entry : entries) {
            std::visit(
                overloaded{
                    [this, &basis, &sort_by_quantum_number_m, &op1, &op2,
                     &delta](const typename GreenTensor<Scalar>::ConstantEntry &ce) {
                        this->hamiltonian->get_matrix() += ce.val() *
                            utils::calculate_tensor_product(basis, basis, op1[ce.row()],
                                                            op2[ce.col()]);
                        if (ce.row() != ce.col() + delta) {
                            sort_by_quantum_number_m = false;
                        }
                    },

                    [this, &basis, &energies, &sort_by_quantum_number_m, &op1, &op2,
                     &delta](const typename GreenTensor<Scalar>::OmegaDependentEntry &oe) {
                        auto tensor_product = utils::calculate_tensor_product(
                            basis, basis, op1[oe.row()], op2[oe.col()]);
                        for (int k = 0; k < tensor_product.outerSize(); ++k) {
                            for (typename Eigen::SparseMatrix<
                                     Scalar, Eigen::RowMajor>::InnerIterator it(tensor_product, k);
                                 it; ++it) {
                                it.valueRef() *= oe.val(energies(it.row()) - energies(it.col()));
                            }
                        }
                        this->hamiltonian->get_matrix() += tensor_product;
                        if (oe.row() != oe.col() + delta) {
                            sort_by_quantum_number_m = false;
                        }
                    }},
                entry);

            this->hamiltonian_is_diagonal = false;
            sort_by_quantum_number_f = false;
        }
    };

    // Dipole-dipole interaction
    add_interaction(green_tensor_ptr->get_entries(1, 1), op.d1, op.d2, 0);

    // Dipole-quadrupole interaction
    add_interaction(green_tensor_ptr->get_entries(1, 2), op.d1, op.q2, -1);

    // Quadrupole-dipole interaction
    add_interaction(green_tensor_ptr->get_entries(2, 1), op.q1, op.d2, +1);

    // Quadrupole-quadrupole interaction
    add_interaction(green_tensor_ptr->get_entries(2, 2), op.q1, op.q2, 0);

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
