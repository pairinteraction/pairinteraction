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
SystemCombined<Scalar> &SystemCombined<Scalar>::set_distance(real_t distance) {
    return set_distance_vector({0, 0, distance});
}

template <typename Scalar>
SystemCombined<Scalar> &
SystemCombined<Scalar>::set_distance_vector(const std::array<real_t, 3> &vector) {
    // https://doi.org/10.1103/PhysRevA.96.062509
    // https://en.wikipedia.org/wiki/Table_of_spherical_harmonics

    this->hamiltonian_requires_construction = true;

    real_t numerical_precision = 10 * std::numeric_limits<real_t>::epsilon();

    if (!traits::NumTraits<Scalar>::is_complex_v && std::abs(vector[1]) > numerical_precision) {
        throw std::invalid_argument(
            "The distance vector must not have a y-component if the scalar type is real.");
    }

    auto to_spherical1 =
        pairinteraction::spherical::CARTESIAN_TO_SPHERICAL1.template cast<std::complex<real_t>>();
    auto to_spherical2 =
        pairinteraction::spherical::CARTESIAN_TO_SPHERICAL2.template cast<std::complex<real_t>>();
    Eigen::Map<const Eigen::Vector3<real_t>> vector_map(vector.data(), vector.size());
    real_t distance = vector_map.norm();
    Eigen::Vector3<real_t> vector_normalized = vector_map / distance;

    SPDLOG_DEBUG("Interatomic distance: {}", distance);

    // If the distance is zero, disable interaction
    if (distance < numerical_precision) {
        green_function_dipole_dipole.setZero();
        green_function_dipole_quadrupole.setZero();
        green_function_quadrupole_dipole.setZero();
        green_function_quadrupole_quadrupole.setZero();
        return *this;
    }

    // Dyadic green function of dipole-dipole interaction
    {
        Eigen::Matrix3<real_t> green_function_cartesian = Eigen::Matrix3<real_t>::Identity() -
            3 * vector_normalized * vector_normalized.transpose();

        auto tmp = to_spherical1 * green_function_cartesian * to_spherical1.adjoint();
        if constexpr (traits::NumTraits<Scalar>::is_complex_v) {
            green_function_dipole_dipole =
                tmp.sparseView(numerical_precision, 1) / std::pow(distance, 3);
        } else {
            green_function_dipole_dipole =
                tmp.real().sparseView(numerical_precision, 1) / std::pow(distance, 3);
            assert(tmp.imag().norm() < numerical_precision);
        }
    }

    SPDLOG_DEBUG("Green function of dipole-dipole interaction:\n{}",
                 fmt::streamed(green_function_dipole_dipole * std::pow(distance, 3)));

    // Dyadic green function of dipole-quadrupole interaction
    {
        Eigen::Matrix<real_t, 3, 9> green_function_cartesian;
        green_function_cartesian.setZero(); // TODO

        auto tmp = to_spherical1 * green_function_cartesian * to_spherical2.adjoint();
        if constexpr (traits::NumTraits<Scalar>::is_complex_v) {
            green_function_dipole_quadrupole =
                tmp.sparseView(numerical_precision, 1) / std::pow(distance, 4);
        } else {
            green_function_dipole_quadrupole =
                tmp.real().sparseView(numerical_precision, 1) / std::pow(distance, 4);
            assert(tmp.imag().norm() < numerical_precision);
        }
    }

    SPDLOG_DEBUG("Green function of dipole-quadrupole interaction:\n{}",
                 fmt::streamed(green_function_dipole_quadrupole * std::pow(distance, 4)));

    // Dyadic green function of quadrupole-dipole interaction
    {
        Eigen::Matrix<real_t, 9, 3> green_function_cartesian;
        green_function_cartesian.setZero(); // TODO

        auto tmp = to_spherical2 * green_function_cartesian * to_spherical1.adjoint();
        if constexpr (traits::NumTraits<Scalar>::is_complex_v) {
            green_function_quadrupole_dipole =
                tmp.sparseView(numerical_precision, 1) / std::pow(distance, 4);
        } else {
            green_function_quadrupole_dipole =
                tmp.real().sparseView(numerical_precision, 1) / std::pow(distance, 4);
            assert(tmp.imag().norm() < numerical_precision);
        }
    }

    SPDLOG_DEBUG("Green function of quadrupole-dipole interaction:\n{}",
                 fmt::streamed(green_function_quadrupole_dipole * std::pow(distance, 4)));

    // Dyadic green function of quadrupole-quadrupole interaction
    {
        Eigen::Matrix<real_t, 9, 9> green_function_cartesian;
        green_function_cartesian.setZero(); // TODO

        auto tmp = to_spherical2 * green_function_cartesian * to_spherical2.adjoint();
        if constexpr (traits::NumTraits<Scalar>::is_complex_v) {
            green_function_quadrupole_quadrupole =
                tmp.sparseView(numerical_precision, 1) / std::pow(distance, 5);
        } else {
            green_function_quadrupole_quadrupole =
                tmp.real().sparseView(numerical_precision, 1) / std::pow(distance, 5);
            assert(tmp.imag().norm() < numerical_precision);
        }
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
            for (typename Eigen::SparseMatrix<Scalar, Eigen::RowMajor>::InnerIterator it1(
                     green_function_dipole_dipole, row);
                 it1; ++it1) {
                this->hamiltonian->get_matrix() +=
                    it1.value() * calculate_tensor_product(basis, d1[it1.row()], d2[it1.col()]);
                if (it1.row() != it1.col()) {
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
            for (typename Eigen::SparseMatrix<Scalar, Eigen::RowMajor>::InnerIterator it1(
                     green_function_dipole_quadrupole, row);
                 it1; ++it1) {
                this->hamiltonian->get_matrix() +=
                    it1.value() * calculate_tensor_product(basis, d1[it1.row()], q2[it1.col()]);
                if (it1.row() != it1.col() - 1) {
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
            for (typename Eigen::SparseMatrix<Scalar, Eigen::RowMajor>::InnerIterator it1(
                     green_function_quadrupole_dipole, row);
                 it1; ++it1) {
                this->hamiltonian->get_matrix() +=
                    it1.value() * calculate_tensor_product(basis, q1[it1.row()], d2[it1.col()]);
                if (it1.row() - 1 != it1.col()) {
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
            for (typename Eigen::SparseMatrix<Scalar, Eigen::RowMajor>::InnerIterator it1(
                     green_function_quadrupole_quadrupole, row);
                 it1; ++it1) {
                this->hamiltonian->get_matrix() +=
                    it1.value() * calculate_tensor_product(basis, q1[it1.row()], q2[it1.col()]);
                if (it1.row() != it1.col()) {
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
                            size_t col_ket_id = col1 * number_of_states2 + col2;
                            if (!basis->has_ket_index(col_ket_id)) {
                                continue;
                            }
                            if (col2 >= static_cast<Eigen::Index>(range_col2.max())) {
                                break;
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
