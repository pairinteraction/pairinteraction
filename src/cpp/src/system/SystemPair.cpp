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
struct GreenFunctions {
    Eigen::SparseMatrix<Scalar, Eigen::RowMajor> dipole_dipole{3, 3};
    Eigen::SparseMatrix<Scalar, Eigen::RowMajor> dipole_quadrupole{3, 6};
    Eigen::SparseMatrix<Scalar, Eigen::RowMajor> quadrupole_dipole{6, 3};
    Eigen::SparseMatrix<Scalar, Eigen::RowMajor> quadrupole_quadrupole{6, 6};
};

template <typename Scalar>
struct OperatorMatrices {
    std::vector<Eigen::SparseMatrix<Scalar, Eigen::RowMajor>> d1;
    std::vector<Eigen::SparseMatrix<Scalar, Eigen::RowMajor>> d2;
    std::vector<Eigen::SparseMatrix<Scalar, Eigen::RowMajor>> q1;
    std::vector<Eigen::SparseMatrix<Scalar, Eigen::RowMajor>> q2;
};

template <typename Scalar>
GreenFunctions<Scalar> construct_green_functions(
    const std::array<typename traits::NumTraits<Scalar>::real_t, 3> &distance_vector, int order) {
    // https://doi.org/10.1103/PhysRevA.96.062509
    // https://doi.org/10.1103/PhysRevA.82.010901
    // https://en.wikipedia.org/wiki/Table_of_spherical_harmonics

    using real_t = typename traits::NumTraits<Scalar>::real_t;
    constexpr real_t numerical_precision = 100 * std::numeric_limits<real_t>::epsilon();

    GreenFunctions<Scalar> green_functions;

    // Normalize the distance vector, return zero green functions if the distance is infinity
    Eigen::Map<const Eigen::Vector3<real_t>> vector_map(distance_vector.data(),
                                                        distance_vector.size());
    real_t distance = vector_map.norm();
    SPDLOG_DEBUG("Interatomic distance: {}", distance);
    if (distance == std::numeric_limits<real_t>::infinity()) {
        return green_functions;
    }
    if (distance < numerical_precision) {
        throw std::invalid_argument("The distance must be greater than zero.");
    }
    if (!traits::NumTraits<Scalar>::is_complex_v &&
        std::abs(distance_vector[1]) > numerical_precision) {
        throw std::invalid_argument(
            "The distance vector must not have a y-component if the scalar type is real.");
    }
    Eigen::Vector3<real_t> vector_normalized = vector_map / distance;

    // Conversion matrices from Cartesian to spherical coordinates
    const auto &to_spherical_kappa1 = pairinteraction::spherical::CARTESIAN_TO_SPHERICAL_KAPPA1
                                          .template cast<std::complex<real_t>>();
    const auto &to_spherical_kappa2 = pairinteraction::spherical::CARTESIAN_TO_SPHERICAL_KAPPA2
                                          .template cast<std::complex<real_t>>();

    // Dyadic green function of dipole-dipole interaction
    if (order >= 3) {
        Eigen::Matrix3<real_t> green_function_cartesian = Eigen::Matrix3<real_t>::Identity() -
            3 * vector_normalized * vector_normalized.transpose();

        auto tmp = to_spherical_kappa1 * green_function_cartesian * to_spherical_kappa1.adjoint();
        if constexpr (traits::NumTraits<Scalar>::is_complex_v) {
            green_functions.dipole_dipole =
                tmp.sparseView(numerical_precision, 1) / std::pow(distance, 3);
        } else {
            green_functions.dipole_dipole =
                tmp.real().sparseView(numerical_precision, 1) / std::pow(distance, 3);
            assert(tmp.imag().norm() < numerical_precision);
        }
    } else {
        green_functions.dipole_dipole.setZero();
    }

    SPDLOG_DEBUG("Green function of dipole-dipole interaction:\n{}",
                 fmt::streamed(green_functions.dipole_dipole * std::pow(distance, 3)));

    // Dyadic green function of dipole-quadrupole interaction
    if (order >= 4) {
        Eigen::Matrix<real_t, 3, 9> green_function_cartesian = Eigen::Matrix<real_t, 3, 9>::Zero();
        for (Eigen::Index q = 0; q < 3; ++q) {
            Eigen::Index row = q;
            for (Eigen::Index j = 0; j < 3; ++j) {
                for (Eigen::Index i = 0; i < 3; ++i) {
                    Eigen::Index col = 3 * j + i;
                    green_function_cartesian(row, col) +=
                        15 * vector_normalized[q] * vector_normalized[j] * vector_normalized[i];
                    if (i == j) {
                        green_function_cartesian(row, col) += -3 * vector_normalized[q];
                    }
                    if (i == q) {
                        green_function_cartesian(row, col) += -3 * vector_normalized[j];
                    }
                    if (j == q) {
                        green_function_cartesian(row, col) += -3 * vector_normalized[i];
                    }
                }
            }
        }

        auto tmp = to_spherical_kappa1 * green_function_cartesian * to_spherical_kappa2.adjoint();
        if constexpr (traits::NumTraits<Scalar>::is_complex_v) {
            green_functions.dipole_quadrupole =
                tmp.sparseView(numerical_precision, 1) / std::pow(distance, 4);
        } else {
            green_functions.dipole_quadrupole =
                tmp.real().sparseView(numerical_precision, 1) / std::pow(distance, 4);
            assert(tmp.imag().norm() < numerical_precision);
        }
    } else {
        green_functions.dipole_quadrupole.setZero();
    }

    SPDLOG_DEBUG("Green function of dipole-quadrupole interaction:\n{}",
                 fmt::streamed(green_functions.dipole_quadrupole * std::pow(distance, 4)));

    // Dyadic green function of quadrupole-dipole interaction
    if (order >= 4) {
        Eigen::Matrix<real_t, 9, 3> green_function_cartesian = Eigen::Matrix<real_t, 9, 3>::Zero();
        for (Eigen::Index q = 0; q < 3; ++q) {
            for (Eigen::Index j = 0; j < 3; ++j) {
                Eigen::Index row = 3 * q + j;
                for (Eigen::Index i = 0; i < 3; ++i) {
                    Eigen::Index col = i;
                    green_function_cartesian(row, col) +=
                        -15 * vector_normalized[q] * vector_normalized[j] * vector_normalized[i];
                    if (i == j) {
                        green_function_cartesian(row, col) += 3 * vector_normalized[q];
                    }
                    if (i == q) {
                        green_function_cartesian(row, col) += 3 * vector_normalized[j];
                    }
                    if (j == q) {
                        green_function_cartesian(row, col) += 3 * vector_normalized[i];
                    }
                }
            }
        }

        auto tmp = to_spherical_kappa2 * green_function_cartesian * to_spherical_kappa1.adjoint();
        if constexpr (traits::NumTraits<Scalar>::is_complex_v) {
            green_functions.quadrupole_dipole =
                tmp.sparseView(numerical_precision, 1) / std::pow(distance, 4);
        } else {
            green_functions.quadrupole_dipole =
                tmp.real().sparseView(numerical_precision, 1) / std::pow(distance, 4);
            assert(tmp.imag().norm() < numerical_precision);
        }
    } else {
        green_functions.quadrupole_dipole.setZero();
    }

    SPDLOG_DEBUG("Green function of quadrupole-dipole interaction:\n{}",
                 fmt::streamed(green_functions.quadrupole_dipole * std::pow(distance, 4)));

    // Dyadic green function of quadrupole-quadrupole interaction
    if (order >= 5) {
        Eigen::Matrix<real_t, 9, 9> green_function_cartesian = Eigen::Matrix<real_t, 9, 9>::Zero();
        for (Eigen::Index q = 0; q < 3; ++q) {
            for (Eigen::Index j = 0; j < 3; ++j) {
                Eigen::Index row = 3 * q + j;
                for (Eigen::Index i = 0; i < 3; ++i) {
                    for (Eigen::Index k = 0; k < 3; ++k) {
                        Eigen::Index col = 3 * i + k;
                        green_function_cartesian(row, col) += 105 * vector_normalized[q] *
                            vector_normalized[j] * vector_normalized[i] * vector_normalized[k];
                        if (i == j) {
                            green_function_cartesian(row, col) +=
                                -15 * vector_normalized[q] * vector_normalized[k];
                        }
                        if (i == q) {
                            green_function_cartesian(row, col) +=
                                -15 * vector_normalized[j] * vector_normalized[k];
                        }
                        if (j == q) {
                            green_function_cartesian(row, col) +=
                                -15 * vector_normalized[i] * vector_normalized[k];
                        }
                        if (k == q) {
                            green_function_cartesian(row, col) +=
                                -15 * vector_normalized[j] * vector_normalized[i];
                        }
                        if (k == j) {
                            green_function_cartesian(row, col) +=
                                -15 * vector_normalized[q] * vector_normalized[i];
                        }
                        if (k == i) {
                            green_function_cartesian(row, col) +=
                                -15 * vector_normalized[q] * vector_normalized[j];
                        }
                        if (q == k && i == j) {
                            green_function_cartesian(row, col) += 3;
                        }
                        if (i == k && j == q) {
                            green_function_cartesian(row, col) += 3;
                        }
                        if (j == k && i == q) {
                            green_function_cartesian(row, col) += 3;
                        }
                    }
                }
            }
        }

        auto tmp = to_spherical_kappa2 * green_function_cartesian * to_spherical_kappa2.adjoint();
        if constexpr (traits::NumTraits<Scalar>::is_complex_v) {
            green_functions.quadrupole_quadrupole =
                tmp.sparseView(numerical_precision, 1) / std::pow(distance, 5);
        } else {
            green_functions.quadrupole_quadrupole =
                tmp.real().sparseView(numerical_precision, 1) / std::pow(distance, 5);
            assert(tmp.imag().norm() < numerical_precision);
        }
    } else {
        green_functions.quadrupole_quadrupole.setZero();
    }

    SPDLOG_DEBUG("Green function of quadrupole-quadrupole interaction:\n{}",
                 fmt::streamed(green_functions.quadrupole_quadrupole * std::pow(distance, 5)));

    return green_functions;
}

template <typename Scalar>
OperatorMatrices<Scalar>
construct_operator_matrices(const GreenFunctions<Scalar> &green_functions,
                            const std::shared_ptr<const BasisAtom<Scalar>> &basis1,
                            const std::shared_ptr<const BasisAtom<Scalar>> &basis2) {
    OperatorMatrices<Scalar> op;

    if (green_functions.dipole_dipole.nonZeros() > 0 ||
        green_functions.dipole_quadrupole.nonZeros() > 0) {
        op.d1.push_back(
            -OperatorAtom<Scalar>(basis1, OperatorType::ELECTRIC_DIPOLE, 1).get_matrix());
        op.d1.push_back(
            OperatorAtom<Scalar>(basis1, OperatorType::ELECTRIC_DIPOLE, 0).get_matrix());
        op.d1.push_back(
            -OperatorAtom<Scalar>(basis1, OperatorType::ELECTRIC_DIPOLE, -1).get_matrix());
    }

    if (green_functions.dipole_dipole.nonZeros() > 0 ||
        green_functions.quadrupole_dipole.nonZeros() > 0) {
        op.d2.push_back(
            OperatorAtom<Scalar>(basis2, OperatorType::ELECTRIC_DIPOLE, -1).get_matrix());
        op.d2.push_back(
            OperatorAtom<Scalar>(basis2, OperatorType::ELECTRIC_DIPOLE, 0).get_matrix());
        op.d2.push_back(
            OperatorAtom<Scalar>(basis2, OperatorType::ELECTRIC_DIPOLE, 1).get_matrix());
    }

    if (green_functions.quadrupole_quadrupole.nonZeros() > 0 ||
        green_functions.quadrupole_dipole.nonZeros() > 0) {
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

    if (green_functions.quadrupole_quadrupole.nonZeros() > 0 ||
        green_functions.dipole_quadrupole.nonZeros() > 0) {
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
SystemPair<Scalar> &SystemPair<Scalar>::set_distance(real_t value) {
    return set_distance_vector({0, 0, value});
}

template <typename Scalar>
SystemPair<Scalar> &SystemPair<Scalar>::set_distance_vector(const std::array<real_t, 3> &vector) {
    this->hamiltonian_requires_construction = true;
    distance_vector = vector;
    return *this;
}

template <typename Scalar>
void SystemPair<Scalar>::construct_hamiltonian() const {
    auto basis = this->hamiltonian->get_basis();
    auto basis1 = basis->get_basis1();
    auto basis2 = basis->get_basis2();

    auto green_functions = construct_green_functions<Scalar>(distance_vector, order);
    auto op = construct_operator_matrices(green_functions, basis1, basis2);

    // Construct the unperturbed Hamiltonian
    this->hamiltonian = std::make_unique<OperatorPair<Scalar>>(basis, OperatorType::ENERGY);
    this->hamiltonian_is_diagonal = true;
    bool sort_by_quantum_number_f = basis->has_quantum_number_f();
    bool sort_by_quantum_number_m = basis->has_quantum_number_m();
    bool sort_by_parity = basis->has_parity();

    // Dipole-dipole interaction
    if (green_functions.dipole_dipole.nonZeros() > 0) {
        for (Eigen::Index row = 0; row < green_functions.dipole_dipole.rows(); ++row) {
            for (typename Eigen::SparseMatrix<Scalar, Eigen::RowMajor>::InnerIterator it(
                     green_functions.dipole_dipole, row);
                 it; ++it) {
                this->hamiltonian->get_matrix() += it.value() *
                    utils::calculate_tensor_product(basis, basis, op.d1[it.row()], op.d2[it.col()]);
                if (it.row() != it.col()) {
                    sort_by_quantum_number_m = false;
                }
            }
        }
        this->hamiltonian_is_diagonal = false;
        sort_by_quantum_number_f = false;
    }

    // Dipole-quadrupole interaction
    if (green_functions.dipole_quadrupole.nonZeros() > 0) {
        for (Eigen::Index row = 0; row < green_functions.dipole_quadrupole.rows(); ++row) {
            for (typename Eigen::SparseMatrix<Scalar, Eigen::RowMajor>::InnerIterator it(
                     green_functions.dipole_quadrupole, row);
                 it; ++it) {
                this->hamiltonian->get_matrix() += it.value() *
                    utils::calculate_tensor_product(basis, basis, op.d1[it.row()], op.q2[it.col()]);
                if (it.row() != it.col() - 1) {
                    sort_by_quantum_number_m = false;
                }
            }
        }
        this->hamiltonian_is_diagonal = false;
        sort_by_quantum_number_f = false;
    }

    // Quadrupole-dipole interaction
    if (green_functions.quadrupole_dipole.nonZeros() > 0) {
        for (Eigen::Index row = 0; row < green_functions.quadrupole_dipole.rows(); ++row) {
            for (typename Eigen::SparseMatrix<Scalar, Eigen::RowMajor>::InnerIterator it(
                     green_functions.quadrupole_dipole, row);
                 it; ++it) {
                this->hamiltonian->get_matrix() += it.value() *
                    utils::calculate_tensor_product(basis, basis, op.q1[it.row()], op.d2[it.col()]);
                if (it.row() - 1 != it.col()) {
                    sort_by_quantum_number_m = false;
                }
            }
        }
        this->hamiltonian_is_diagonal = false;
        sort_by_quantum_number_f = false;
    }

    // Quadrupole-quadrupole interaction
    if (green_functions.quadrupole_quadrupole.nonZeros() > 0) {
        for (Eigen::Index row = 0; row < green_functions.quadrupole_quadrupole.rows(); ++row) {
            for (typename Eigen::SparseMatrix<Scalar, Eigen::RowMajor>::InnerIterator it(
                     green_functions.quadrupole_quadrupole, row);
                 it; ++it) {
                this->hamiltonian->get_matrix() += it.value() *
                    utils::calculate_tensor_product(basis, basis, op.q1[it.row()], op.q2[it.col()]);
                if (it.row() != it.col()) {
                    sort_by_quantum_number_m = false;
                }
            }
        }
        this->hamiltonian_is_diagonal = false;
        sort_by_quantum_number_f = false;
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
