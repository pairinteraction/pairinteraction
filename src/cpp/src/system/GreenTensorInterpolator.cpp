// SPDX-FileCopyrightText: 2025 PairInteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "pairinteraction/system/GreenTensorInterpolator.hpp"

#include "pairinteraction/utils/eigen_assertion.hpp"
#include "pairinteraction/utils/eigen_compat.hpp"
#include "pairinteraction/utils/spherical.hpp"
#include "pairinteraction/utils/traits.hpp"

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <array>
#include <cmath>
#include <complex>
#include <map>
#include <spdlog/spdlog.h>
#include <unsupported/Eigen/Splines>

namespace pairinteraction {
template <typename Scalar>
GreenTensorInterpolator<Scalar>::ConstantEntry::ConstantEntry(int row, int col, Scalar val)
    : row_(row), col_(col), val_(val) {}

template <typename Scalar>
Scalar GreenTensorInterpolator<Scalar>::ConstantEntry::val() const {
    return val_;
}

template <typename Scalar>
int GreenTensorInterpolator<Scalar>::ConstantEntry::row() const noexcept {
    return row_;
}

template <typename Scalar>
int GreenTensorInterpolator<Scalar>::ConstantEntry::col() const noexcept {
    return col_;
}

template <typename Scalar>
GreenTensorInterpolator<Scalar>::OmegaDependentEntry::OmegaDependentEntry(
    int row, int col, Eigen::Spline<real_t, 1> real_spline, Eigen::Spline<real_t, 1> imag_spline)
    : row_(row), col_(col), real_spline(std::move(real_spline)),
      imag_spline(std::move(imag_spline)) {}

template <typename Scalar>
Scalar GreenTensorInterpolator<Scalar>::OmegaDependentEntry::val(double omega) const {
    if constexpr (traits::NumTraits<Scalar>::is_complex_v) {
        return {real_spline(omega)(0), imag_spline(omega)(0)};
    } else {
        return real_spline(omega)(0);
    }
}

template <typename Scalar>
int GreenTensorInterpolator<Scalar>::OmegaDependentEntry::row() const noexcept {
    return row_;
}

template <typename Scalar>
int GreenTensorInterpolator<Scalar>::OmegaDependentEntry::col() const noexcept {
    return col_;
}

template <typename Scalar>
void GreenTensorInterpolator<Scalar>::create_entries_from_cartesian(
    int kappa1, int kappa2, const Eigen::MatrixX<Scalar> &tensor_in_cartesian_coordinates) {

    const real_t scale = tensor_in_cartesian_coordinates.norm();
    const real_t numerical_precision = 100 * scale * std::numeric_limits<real_t>::epsilon();

    Eigen::SparseMatrix<complex_t> tensor =
        (spherical::get_transformator<complex_t>(kappa1) * tensor_in_cartesian_coordinates *
         spherical::get_transformator<complex_t>(kappa2).adjoint())
            .sparseView(1, numerical_precision);

    std::vector<Entry> entries;
    for (int k = 0; k < tensor.outerSize(); ++k) {
        for (typename Eigen::SparseMatrix<complex_t>::InnerIterator it(tensor, k); it; ++it) {
            if constexpr (traits::NumTraits<Scalar>::is_complex_v) {
                entries.emplace_back(ConstantEntry(it.row(), it.col(), it.value()));
            } else {
                entries.emplace_back(ConstantEntry(it.row(), it.col(), it.value().real()));
                assert(abs(it.value().imag()) < numerical_precision);
            }
        }
    }
    entries_map[{kappa1, kappa2}] = std::move(entries);
}

template <typename Scalar>
void GreenTensorInterpolator<Scalar>::create_entries_from_cartesian(
    int kappa1, int kappa2,
    const std::vector<Eigen::MatrixX<Scalar>> &tensors_in_cartesian_coordinates,
    const std::vector<double> &omegas) {

    if (tensors_in_cartesian_coordinates.size() != omegas.size()) {
        throw std::invalid_argument("The number of tensors and omegas must match.");
    }

    if (tensors_in_cartesian_coordinates.size() < 4) {
        throw std::invalid_argument(
            "At least 4 tensors are required for the applied cubic spline interpolation.");
    }

    auto num_knots = static_cast<int>(omegas.size());
    Eigen::Map<const Eigen::RowVectorXd> knots(omegas.data(), num_knots);

    constexpr int spline_degree = 3; // cubic spline interpolation

    // Temporary storage wih key = (row, col) and value = vector of one double per omega
    std::map<std::pair<int, int>, std::pair<Eigen::RowVectorXd, Eigen::RowVectorXd>> temp_map;
    for (int idx = 0; idx < num_knots; ++idx) {

        const real_t scale = tensors_in_cartesian_coordinates[idx].norm();
        const real_t numerical_precision = 100 * scale * std::numeric_limits<real_t>::epsilon();

        Eigen::SparseMatrix<complex_t> tensor =
            (spherical::get_transformator<complex_t>(kappa1) *
             tensors_in_cartesian_coordinates[idx] *
             spherical::get_transformator<complex_t>(kappa2).adjoint())
                .sparseView(1, numerical_precision);

        for (int k = 0; k < tensor.outerSize(); ++k) {
            for (typename Eigen::SparseMatrix<complex_t>::InnerIterator it(tensor, k); it; ++it) {
                std::pair<int, int> key{it.row(), it.col()};
                auto &[vec_real, vec_imag] =
                    temp_map
                        .try_emplace(key, Eigen::RowVectorXd::Zero(num_knots),
                                     Eigen::RowVectorXd::Zero(num_knots))
                        .first->second;
                vec_real(idx) = it.value().real();
                if constexpr (traits::NumTraits<Scalar>::is_complex_v) {
                    vec_imag(idx) = it.value().imag();
                } else {
                    assert(abs(it.value().imag()) < numerical_precision);
                }
            }
        }
    }

    // Set the green tensor entries with spline interpolation
    std::vector<Entry> entries;
    entries.reserve(temp_map.size());
    for (const auto &[key, value] : temp_map) {
        const auto &[vec_real, vec_imag] = value;
        const auto &[row, col] = key;

        Eigen::Spline<real_t, 1> real_spline =
            Eigen::SplineFitting<Eigen::Spline<real_t, 1>>::Interpolate(vec_real, spline_degree,
                                                                        knots);

        Eigen::Spline<real_t, 1> imag_spline;
        if constexpr (traits::NumTraits<Scalar>::is_complex_v) {
            imag_spline = Eigen::SplineFitting<Eigen::Spline<real_t, 1>>::Interpolate(
                vec_imag, spline_degree, knots);
        }

        entries.emplace_back(
            OmegaDependentEntry(row, col, std::move(real_spline), std::move(imag_spline)));
    }
    entries_map[{kappa1, kappa2}] = std::move(entries);
}

template <typename Scalar>
const std::vector<typename GreenTensorInterpolator<Scalar>::Entry> &
GreenTensorInterpolator<Scalar>::get_spherical_entries(int kappa1, int kappa2) const {
    if (auto it = entries_map.find({kappa1, kappa2}); it != entries_map.end()) {
        return it->second;
    }
    static const std::vector<Entry> empty_entries;
    return empty_entries;
}

template <typename Scalar>
GreenTensorInterpolator<Scalar> GreenTensorInterpolator<Scalar>::from_multipole_expansion(
    const std::array<real_t, 3> &distance_vector, int interaction_order, int minimal_kappa1,
    int minimal_kappa2) {
    // https://doi.org/10.1103/PhysRevA.96.062509
    // https://doi.org/10.1103/PhysRevA.82.010901
    // https://en.wikipedia.org/wiki/Table_of_spherical_harmonics

    GreenTensorInterpolator<Scalar> green_tensor_interpolator;

    // Normalize the distance vector, return zero green tensor if the distance is infinity
    Eigen::Map<const Eigen::Vector3<real_t>> vector_map(distance_vector.data(),
                                                        distance_vector.size());
    real_t distance = vector_map.norm();
    SPDLOG_DEBUG("Interatomic distance: {}", distance);
    if (!std::isfinite(distance)) {
        return green_tensor_interpolator;
    }
    Eigen::Vector3<real_t> unitvec = vector_map / distance;

    auto is_required = [interaction_order, minimal_kappa1, minimal_kappa2](int kappa1, int kappa2) {
        return kappa1 + kappa2 + 1 <= interaction_order && kappa1 >= minimal_kappa1 &&
            kappa2 >= minimal_kappa2;
    };

    // The cartesian tensors follow from the Taylor expansion of the Coulomb interaction
    // 1/|R + r2 - r1| in the coordinates r1 and r2 of the two charge distributions. The tensor
    // for the orders (kappa1, kappa2) is (-1)^kappa1 * d^(kappa1+kappa2)/dR^(kappa1+kappa2) (1/R).
    // Note that the factor 1/(kappa1! * kappa2!) of the Taylor expansion is not included because
    // it is absorbed in the normalization of the cartesian-to-spherical transformation, see
    // spherical.hpp.

    // Green function of monopole-monopole interaction
    if (is_required(0, 0)) {
        Eigen::Matrix<real_t, 1, 1> entries = Eigen::Matrix<real_t, 1, 1>::Constant(1);

        green_tensor_interpolator.create_entries_from_cartesian(
            0, 0, (entries / distance).template cast<Scalar>());
    }

    // Green function of monopole-dipole interaction
    if (is_required(0, 1)) {
        Eigen::Matrix<real_t, 1, 3> entries = -unitvec.transpose();

        green_tensor_interpolator.create_entries_from_cartesian(
            0, 1, (entries / std::pow(distance, 2)).template cast<Scalar>());
    }

    // Green function of dipole-monopole interaction
    if (is_required(1, 0)) {
        const Eigen::Matrix<real_t, 3, 1> &entries = unitvec;

        green_tensor_interpolator.create_entries_from_cartesian(
            1, 0, (entries / std::pow(distance, 2)).template cast<Scalar>());
    }

    // Green function of monopole-quadrupole interaction
    if (is_required(0, 2)) {
        Eigen::Matrix<real_t, 1, 9> entries = Eigen::Matrix<real_t, 1, 9>::Zero();
        for (Eigen::Index j = 0; j < 3; ++j) {
            for (Eigen::Index i = 0; i < 3; ++i) {
                Eigen::Index col = 3 * j + i;
                real_t v = 3 * unitvec[j] * unitvec[i];
                if (i == j) v += -1;
                entries(0, col) += v;
            }
        }

        green_tensor_interpolator.create_entries_from_cartesian(
            0, 2, (entries / std::pow(distance, 3)).template cast<Scalar>());
    }

    // Green function of quadrupole-monopole interaction
    if (is_required(2, 0)) {
        Eigen::Matrix<real_t, 9, 1> entries = Eigen::Matrix<real_t, 9, 1>::Zero();
        for (Eigen::Index q = 0; q < 3; ++q) {
            for (Eigen::Index j = 0; j < 3; ++j) {
                Eigen::Index row = 3 * q + j;
                real_t v = 3 * unitvec[q] * unitvec[j];
                if (q == j) v += -1;
                entries(row, 0) += v;
            }
        }

        green_tensor_interpolator.create_entries_from_cartesian(
            2, 0, (entries / std::pow(distance, 3)).template cast<Scalar>());
    }

    // Dyadic green function of dipole-dipole interaction
    if (is_required(1, 1)) {
        Eigen::Matrix3<Scalar> entries =
            Eigen::Matrix3<real_t>::Identity() - 3 * unitvec * unitvec.transpose();

        green_tensor_interpolator.create_entries_from_cartesian(
            1, 1, (entries / std::pow(distance, 3)).template cast<Scalar>());
    }

    // Dyadic green function of dipole-quadrupole interaction
    if (is_required(1, 2)) {
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

        green_tensor_interpolator.create_entries_from_cartesian(
            1, 2, (entries / std::pow(distance, 4)).template cast<Scalar>());
    }

    // Dyadic green function of quadrupole-dipole interaction
    if (is_required(2, 1)) {
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

        green_tensor_interpolator.create_entries_from_cartesian(
            2, 1, (entries / std::pow(distance, 4)).template cast<Scalar>());
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

        green_tensor_interpolator.create_entries_from_cartesian(
            2, 2, (entries / std::pow(distance, 5)).template cast<Scalar>());
    }

    return green_tensor_interpolator;
}

// Explicit instantiations
template class GreenTensorInterpolator<double>;
template class GreenTensorInterpolator<std::complex<double>>;
} // namespace pairinteraction
