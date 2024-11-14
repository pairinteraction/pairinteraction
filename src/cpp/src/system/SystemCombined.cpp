#include "pairinteraction/system/SystemCombined.hpp"

#include "pairinteraction/basis/BasisAtom.hpp"
#include "pairinteraction/basis/BasisCombined.hpp"
#include "pairinteraction/enums/OperatorType.hpp"
#include "pairinteraction/enums/Parity.hpp"
#include "pairinteraction/enums/TransformationType.hpp"
#include "pairinteraction/ket/KetAtom.hpp"
#include "pairinteraction/ket/KetCombined.hpp"
#include "pairinteraction/operator/OperatorAtom.hpp"
#include "pairinteraction/operator/OperatorCombined.hpp"
#include "pairinteraction/system/SystemAtom.hpp"
#include "pairinteraction/utils/Range.hpp"
#include "pairinteraction/utils/eigen_assertion.hpp"
#include "pairinteraction/utils/eigen_compat.hpp"
#include "pairinteraction/utils/spherical.hpp"
#include "pairinteraction/utils/streamed.hpp"

#include <Eigen/SparseCore>
#include <algorithm>
#include <array>
#include <complex>
#include <limits>
#include <memory>
#include <oneapi/tbb.h>
#include <spdlog/spdlog.h>
#include <vector>

namespace pairinteraction {
template <typename Scalar>
SystemCombined<Scalar>::SystemCombined(std::shared_ptr<const basis_t> basis)
    : System<SystemCombined<Scalar>>(std::move(basis)) {}

template <typename Scalar>
SystemCombined<Scalar> &SystemCombined<Scalar>::set_order(int value) {
    this->hamiltonian_requires_construction = true;
    if (value < 3 || value > 5) {
        throw std::invalid_argument("The order must be 3, 4, or 5.");
    }
    order = value;
    return *this;
}

template <typename Scalar>
SystemCombined<Scalar> &SystemCombined<Scalar>::set_distance(real_t value) {
    return set_distance_vector({0, 0, value});
}

template <typename Scalar>
SystemCombined<Scalar> &
SystemCombined<Scalar>::set_distance_vector(const std::array<real_t, 3> &vector) {
    // https://doi.org/10.1103/PhysRevA.96.062509
    // https://doi.org/10.1103/PhysRevA.82.010901
    // https://en.wikipedia.org/wiki/Table_of_spherical_harmonics

    this->hamiltonian_requires_construction = true;

    real_t numerical_precision = 10 * std::numeric_limits<real_t>::epsilon();

    if (!traits::NumTraits<Scalar>::is_complex_v && std::abs(vector[1]) > numerical_precision) {
        throw std::invalid_argument(
            "The distance vector must not have a y-component if the scalar type is real.");
    }

    const auto &to_spherical_kappa1 = pairinteraction::spherical::CARTESIAN_TO_SPHERICAL_KAPPA1
                                          .template cast<std::complex<real_t>>();
    const auto &to_spherical_kappa2 = pairinteraction::spherical::CARTESIAN_TO_SPHERICAL_KAPPA2
                                          .template cast<std::complex<real_t>>();
    Eigen::Map<const Eigen::Vector3<real_t>> vector_map(vector.data(), vector.size());
    real_t distance = vector_map.norm();
    Eigen::Vector3<real_t> vector_normalized = vector_map / distance;

    SPDLOG_DEBUG("Interatomic distance: {}", distance);

    // Dyadic green function of dipole-dipole interaction
    if (order >= 3) {
        Eigen::Matrix3<real_t> green_function_cartesian = Eigen::Matrix3<real_t>::Identity() -
            3 * vector_normalized * vector_normalized.transpose();

        auto tmp = to_spherical_kappa1 * green_function_cartesian * to_spherical_kappa1.adjoint();
        if constexpr (traits::NumTraits<Scalar>::is_complex_v) {
            green_function_dipole_dipole =
                tmp.sparseView(numerical_precision, 1) / std::pow(distance, 3);
        } else {
            green_function_dipole_dipole =
                tmp.real().sparseView(numerical_precision, 1) / std::pow(distance, 3);
            assert(tmp.imag().norm() < numerical_precision);
        }
    } else {
        green_function_dipole_dipole.setZero();
    }

    SPDLOG_DEBUG("Green function of dipole-dipole interaction:\n{}",
                 fmt::streamed(green_function_dipole_dipole * std::pow(distance, 3)));

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
            green_function_dipole_quadrupole =
                tmp.sparseView(numerical_precision, 1) / std::pow(distance, 4);
        } else {
            green_function_dipole_quadrupole =
                tmp.real().sparseView(numerical_precision, 1) / std::pow(distance, 4);
            assert(tmp.imag().norm() < numerical_precision);
        }
    } else {
        green_function_dipole_quadrupole.setZero();
    }

    SPDLOG_DEBUG("Green function of dipole-quadrupole interaction:\n{}",
                 fmt::streamed(green_function_dipole_quadrupole * std::pow(distance, 4)));

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
            green_function_quadrupole_dipole =
                tmp.sparseView(numerical_precision, 1) / std::pow(distance, 4);
        } else {
            green_function_quadrupole_dipole =
                tmp.real().sparseView(numerical_precision, 1) / std::pow(distance, 4);
            assert(tmp.imag().norm() < numerical_precision);
        }
    } else {
        green_function_quadrupole_dipole.setZero();
    }

    SPDLOG_DEBUG("Green function of quadrupole-dipole interaction:\n{}",
                 fmt::streamed(green_function_quadrupole_dipole * std::pow(distance, 4)));

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
            green_function_quadrupole_quadrupole =
                tmp.sparseView(numerical_precision, 1) / std::pow(distance, 5);
        } else {
            green_function_quadrupole_quadrupole =
                tmp.real().sparseView(numerical_precision, 1) / std::pow(distance, 5);
            assert(tmp.imag().norm() < numerical_precision);
        }
    } else {
        green_function_quadrupole_quadrupole.setZero();
    }

    SPDLOG_DEBUG("Green function of quadrupole-quadrupole interaction:\n{}",
                 fmt::streamed(green_function_quadrupole_quadrupole * std::pow(distance, 5)));

    return *this;
}

template <typename Scalar>
void SystemCombined<Scalar>::construct_hamiltonian() const {
    auto basis = this->hamiltonian->get_basis();
    auto basis1 = basis->get_basis1();
    auto basis2 = basis->get_basis2();

    // Construct the unperturbed Hamiltonian
    this->hamiltonian = std::make_unique<OperatorCombined<Scalar>>(basis, OperatorType::ENERGY);
    this->hamiltonian_is_diagonal = true;
    bool sort_by_quantum_number_f = basis->has_quantum_number_f();
    bool sort_by_quantum_number_m = basis->has_quantum_number_m();
    bool sort_by_parity = basis->has_parity();

    // Get operator matrices
    std::vector<Eigen::SparseMatrix<Scalar, Eigen::RowMajor>> d1;
    std::vector<Eigen::SparseMatrix<Scalar, Eigen::RowMajor>> d2;
    std::vector<Eigen::SparseMatrix<Scalar, Eigen::RowMajor>> q1;
    std::vector<Eigen::SparseMatrix<Scalar, Eigen::RowMajor>> q2;

    if (green_function_dipole_dipole.nonZeros() > 0 ||
        green_function_dipole_quadrupole.nonZeros() > 0) {
        d1.push_back(-OperatorAtom<Scalar>(basis1, OperatorType::ELECTRIC_DIPOLE, 1).get_matrix());
        d1.push_back(OperatorAtom<Scalar>(basis1, OperatorType::ELECTRIC_DIPOLE, 0).get_matrix());
        d1.push_back(-OperatorAtom<Scalar>(basis1, OperatorType::ELECTRIC_DIPOLE, -1).get_matrix());
    }

    if (green_function_dipole_dipole.nonZeros() > 0 ||
        green_function_quadrupole_dipole.nonZeros() > 0) {
        d2.push_back(OperatorAtom<Scalar>(basis2, OperatorType::ELECTRIC_DIPOLE, -1).get_matrix());
        d2.push_back(OperatorAtom<Scalar>(basis2, OperatorType::ELECTRIC_DIPOLE, 0).get_matrix());
        d2.push_back(OperatorAtom<Scalar>(basis2, OperatorType::ELECTRIC_DIPOLE, 1).get_matrix());
    }

    if (green_function_quadrupole_quadrupole.nonZeros() > 0 ||
        green_function_quadrupole_dipole.nonZeros() > 0) {
        q1.push_back(
            OperatorAtom<Scalar>(basis1, OperatorType::ELECTRIC_QUADRUPOLE, 2).get_matrix());
        q1.push_back(
            -OperatorAtom<Scalar>(basis1, OperatorType::ELECTRIC_QUADRUPOLE, 1).get_matrix());
        q1.push_back(
            OperatorAtom<Scalar>(basis1, OperatorType::ELECTRIC_QUADRUPOLE, 0).get_matrix());
        q1.push_back(
            -OperatorAtom<Scalar>(basis1, OperatorType::ELECTRIC_QUADRUPOLE, -1).get_matrix());
        q1.push_back(
            OperatorAtom<Scalar>(basis1, OperatorType::ELECTRIC_QUADRUPOLE, -2).get_matrix());
    }

    if (green_function_quadrupole_quadrupole.nonZeros() > 0 ||
        green_function_dipole_quadrupole.nonZeros() > 0) {
        q2.push_back(
            OperatorAtom<Scalar>(basis2, OperatorType::ELECTRIC_QUADRUPOLE, -2).get_matrix());
        q2.push_back(
            OperatorAtom<Scalar>(basis2, OperatorType::ELECTRIC_QUADRUPOLE, -1).get_matrix());
        q2.push_back(
            OperatorAtom<Scalar>(basis2, OperatorType::ELECTRIC_QUADRUPOLE, 0).get_matrix());
        q2.push_back(
            OperatorAtom<Scalar>(basis2, OperatorType::ELECTRIC_QUADRUPOLE, 1).get_matrix());
        q2.push_back(
            OperatorAtom<Scalar>(basis2, OperatorType::ELECTRIC_QUADRUPOLE, 2).get_matrix());
    }

    // Dipole-dipole interaction
    if (green_function_dipole_dipole.nonZeros() > 0) {
        for (Eigen::Index row = 0; row < green_function_dipole_dipole.rows(); ++row) {
            for (typename Eigen::SparseMatrix<Scalar, Eigen::RowMajor>::InnerIterator it(
                     green_function_dipole_dipole, row);
                 it; ++it) {
                this->hamiltonian->get_matrix() +=
                    it.value() * calculate_tensor_product(basis, d1[it.row()], d2[it.col()]);
                if (it.row() != it.col()) {
                    sort_by_quantum_number_m = false;
                }
            }
        }
        this->hamiltonian_is_diagonal = false;
        sort_by_quantum_number_f = false;
    }

    // Dipole-quadrupole interaction
    if (green_function_dipole_quadrupole.nonZeros() > 0) {
        for (Eigen::Index row = 0; row < green_function_dipole_quadrupole.rows(); ++row) {
            for (typename Eigen::SparseMatrix<Scalar, Eigen::RowMajor>::InnerIterator it(
                     green_function_dipole_quadrupole, row);
                 it; ++it) {
                if (it.col() > 4) {
                    throw std::runtime_error(
                        "p_00 is needed but not yet contained in the database.");
                }
                this->hamiltonian->get_matrix() +=
                    it.value() * calculate_tensor_product(basis, d1[it.row()], q2[it.col()]);
                if (it.row() != it.col() - 1) {
                    sort_by_quantum_number_m = false;
                }
            }
        }
        this->hamiltonian_is_diagonal = false;
        sort_by_quantum_number_f = false;
    }

    // Quadrupole-dipole interaction
    if (green_function_quadrupole_dipole.nonZeros() > 0) {
        for (Eigen::Index row = 0; row < green_function_quadrupole_dipole.rows(); ++row) {
            for (typename Eigen::SparseMatrix<Scalar, Eigen::RowMajor>::InnerIterator it(
                     green_function_quadrupole_dipole, row);
                 it; ++it) {
                if (it.row() > 4) {
                    throw std::runtime_error(
                        "p_00 is needed but not yet contained in the database.");
                }
                this->hamiltonian->get_matrix() +=
                    it.value() * calculate_tensor_product(basis, q1[it.row()], d2[it.col()]);
                if (it.row() - 1 != it.col()) {
                    sort_by_quantum_number_m = false;
                }
            }
        }
        this->hamiltonian_is_diagonal = false;
        sort_by_quantum_number_f = false;
    }

    // Quadrupole-quadrupole interaction
    if (green_function_quadrupole_quadrupole.nonZeros() > 0) {
        for (Eigen::Index row = 0; row < green_function_quadrupole_quadrupole.rows(); ++row) {
            for (typename Eigen::SparseMatrix<Scalar, Eigen::RowMajor>::InnerIterator it(
                     green_function_quadrupole_quadrupole, row);
                 it; ++it) {
                if (it.row() > 4 || it.col() > 4) {
                    throw std::runtime_error(
                        "p_00 is needed but not yet contained in the database.");
                }
                this->hamiltonian->get_matrix() +=
                    it.value() * calculate_tensor_product(basis, q1[it.row()], q2[it.col()]);
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

template <typename Scalar>
Eigen::SparseMatrix<Scalar, Eigen::RowMajor> SystemCombined<Scalar>::calculate_tensor_product(
    const std::shared_ptr<const basis_t> &basis,
    const Eigen::SparseMatrix<Scalar, Eigen::RowMajor> &matrix1,
    const Eigen::SparseMatrix<Scalar, Eigen::RowMajor> &matrix2) {
    real_t numerical_precision = 10 * std::numeric_limits<real_t>::epsilon();

    size_t number_of_states2 = basis->get_basis2()->get_number_of_states();

    oneapi::tbb::concurrent_vector<Eigen::Triplet<Scalar>> triplets;

    // Loop over the rows of the first matrix in parallel (outer index == row)
    oneapi::tbb::parallel_for(
        oneapi::tbb::blocked_range<Eigen::Index>(0, matrix1.outerSize()), [&](const auto &range) {
            for (Eigen::Index row1 = range.begin(); row1 != range.end(); ++row1) {

                const auto &range_row2 = basis->get_index_range(row1);

                // Loop over the rows of the second matrix that are energetically allowed
                for (auto row2 = static_cast<Eigen::Index>(range_row2.min());
                     row2 < static_cast<Eigen::Index>(range_row2.max()); ++row2) {

                    size_t row_ket_id = row1 * number_of_states2 + row2;
                    if (!basis->has_ket_index(row_ket_id)) {
                        continue;
                    }
                    Eigen::Index row = basis->get_ket_index_from_id(row_ket_id);

                    // Loop over the non-zero column elements of the first matrix
                    for (typename Eigen::SparseMatrix<Scalar, Eigen::RowMajor>::InnerIterator it1(
                             matrix1, row1);
                         it1; ++it1) {

                        Eigen::Index col1 = it1.col();
                        Scalar value1 = it1.value();

                        // Calculate the minimum and maximum values of the index pointer of the
                        // second matrix
                        Eigen::Index begin_idxptr2 = matrix2.outerIndexPtr()[row2];
                        Eigen::Index end_idxptr2 = matrix2.outerIndexPtr()[row2 + 1];

                        // The minimum value is chosen such that we start with an energetically
                        // allowed column
                        const auto &range_col2 = basis->get_index_range(it1.index());
                        begin_idxptr2 +=
                            std::distance(matrix2.innerIndexPtr() + begin_idxptr2,
                                          std::lower_bound(matrix2.innerIndexPtr() + begin_idxptr2,
                                                           matrix2.innerIndexPtr() + end_idxptr2,
                                                           range_col2.min()));

                        // Loop over the non-zero column elements of the second matrix that are
                        // energetically allowed (we break the loop if the index pointer corresponds
                        // to a column that is not energetically allowed)
                        for (Eigen::Index idxptr2 = begin_idxptr2; idxptr2 < end_idxptr2;
                             ++idxptr2) {

                            Eigen::Index col2 = matrix2.innerIndexPtr()[idxptr2];
                            if (col2 >= static_cast<Eigen::Index>(range_col2.max())) {
                                break;
                            }
                            size_t col_ket_id = col1 * number_of_states2 + col2;
                            if (!basis->has_ket_index(col_ket_id)) {
                                continue;
                            }
                            Scalar value2 = matrix2.valuePtr()[idxptr2];
                            Eigen::Index col = basis->get_ket_index_from_id(col_ket_id);

                            // Store the entry
                            Scalar value = value1 * value2;
                            if (std::abs(value) > numerical_precision) {
                                triplets.emplace_back(row, col, value);
                            }
                        }
                    }
                }
            }
        });

    // Construct the combined matrix from the triplets
    Eigen::SparseMatrix<Scalar, Eigen::RowMajor> matrix(basis->get_number_of_states(),
                                                        basis->get_number_of_states());
    matrix.setFromTriplets(triplets.begin(), triplets.end());
    matrix.makeCompressed();

    return matrix;
}

// Explicit instantiations
template class SystemCombined<float>;
template class SystemCombined<double>;
template class SystemCombined<std::complex<float>>;
template class SystemCombined<std::complex<double>>;
} // namespace pairinteraction
