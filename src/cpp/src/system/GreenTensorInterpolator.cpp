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

namespace {
// The cartesian entries of the green tensors constructed below are the interaction tensors of
// the multipole expansion of the Coulomb interaction in free space,
// (-1)^kappa1 * R^(kappa1+kappa2+1) * grad^(kappa1+kappa2) (1/R), where the sign accounts for
// the electron coordinate of the first atom entering the interatomic distance with a negative
// sign. Their normalization is fixed such that, together with the cartesian-to-spherical
// transformators of utils/spherical.cpp, the known spherical multipole expansion is
// reproduced (checked in GreenTensorInterpolator.test.cpp). References:
// - S. Weber et al., J. Phys. B 50, 133001 (2017), https://doi.org/10.1088/1361-6455/aa743a
//   [multipole expansion of the Rydberg-Rydberg interaction in spherical harmonics,
//   Eqs. (6)-(8)]
// - A. J. Stone, The Theory of Intermolecular Forces, 2nd ed. (Oxford University Press,
//   2013), https://doi.org/10.1093/acprof:oso/9780199672394.001.0001
//   [explicit cartesian interaction tensors up to rank four, i.e., up to
//   quadrupole-quadrupole and dipole-octupole interaction]
// - J. Block and S. Scheel, Phys. Rev. A 96, 062509 (2017),
//   https://doi.org/10.1103/PhysRevA.96.062509
//   [green tensor formulation of the dipole-dipole interaction]
// - J. A. Crosse et al., Phys. Rev. A 82, 010901(R) (2010),
//   https://doi.org/10.1103/PhysRevA.82.010901
//   [green tensor formulation involving quadrupole transitions]

// Dyadic green function of dipole-dipole interaction
template <typename Real>
Eigen::Matrix3<Real> dipole_dipole_tensor(const Eigen::Vector3<Real> &unitvec) {
    return Eigen::Matrix3<Real>::Identity() - 3 * unitvec * unitvec.transpose();
}

// Dyadic green function of dipole-quadrupole interaction
template <typename Real>
Eigen::Matrix<Real, 3, 9> dipole_quadrupole_tensor(const Eigen::Vector3<Real> &unitvec) {
    Eigen::Matrix<Real, 3, 9> entries = Eigen::Matrix<Real, 3, 9>::Zero();
    for (Eigen::Index q = 0; q < 3; ++q) {
        Eigen::Index row = q;
        for (Eigen::Index j = 0; j < 3; ++j) {
            for (Eigen::Index i = 0; i < 3; ++i) {
                Eigen::Index col = 3 * j + i;
                Real v = 15 * unitvec[q] * unitvec[j] * unitvec[i];
                if (i == j) v += -3 * unitvec[q];
                if (i == q) v += -3 * unitvec[j];
                if (j == q) v += -3 * unitvec[i];
                entries(row, col) += v;
            }
        }
    }
    return entries;
}

// Dyadic green function of quadrupole-dipole interaction
template <typename Real>
Eigen::Matrix<Real, 9, 3> quadrupole_dipole_tensor(const Eigen::Vector3<Real> &unitvec) {
    Eigen::Matrix<Real, 9, 3> entries = Eigen::Matrix<Real, 9, 3>::Zero();
    for (Eigen::Index q = 0; q < 3; ++q) {
        for (Eigen::Index j = 0; j < 3; ++j) {
            Eigen::Index row = 3 * q + j;
            for (Eigen::Index i = 0; i < 3; ++i) {
                Eigen::Index col = i;
                Real v = -15 * unitvec[q] * unitvec[j] * unitvec[i];
                if (i == j) v += 3 * unitvec[q];
                if (i == q) v += 3 * unitvec[j];
                if (j == q) v += 3 * unitvec[i];
                entries(row, col) += v;
            }
        }
    }
    return entries;
}

// Dyadic green function of quadrupole-quadrupole interaction
template <typename Real>
Eigen::Matrix<Real, 9, 9> quadrupole_quadrupole_tensor(const Eigen::Vector3<Real> &unitvec) {
    Eigen::Matrix<Real, 9, 9> entries = Eigen::Matrix<Real, 9, 9>::Zero();
    for (Eigen::Index q = 0; q < 3; ++q) {
        for (Eigen::Index j = 0; j < 3; ++j) {
            Eigen::Index row = 3 * q + j;
            for (Eigen::Index i = 0; i < 3; ++i) {
                for (Eigen::Index k = 0; k < 3; ++k) {
                    Eigen::Index col = 3 * i + k;
                    Real v = 105 * unitvec[q] * unitvec[j] * unitvec[i] * unitvec[k];
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
    return entries;
}

// Common coefficient of the rank-four interaction tensors of the dipole-octupole and
// octupole-dipole interaction. The expression is symmetric under exchanging the roles of the
// indices, so both tensors share it.
template <typename Real>
Real dipole_octupole_coefficient(const Eigen::Vector3<Real> &unitvec, Eigen::Index q,
                                 Eigen::Index j, Eigen::Index i, Eigen::Index k) {
    Real v = -105 * unitvec[q] * unitvec[j] * unitvec[i] * unitvec[k];
    if (i == j) v += 15 * unitvec[q] * unitvec[k];
    if (i == q) v += 15 * unitvec[j] * unitvec[k];
    if (j == q) v += 15 * unitvec[i] * unitvec[k];
    if (k == q) v += 15 * unitvec[j] * unitvec[i];
    if (k == j) v += 15 * unitvec[q] * unitvec[i];
    if (k == i) v += 15 * unitvec[q] * unitvec[j];
    if (q == k && i == j) v += -3;
    if (i == k && j == q) v += -3;
    if (j == k && i == q) v += -3;
    return v;
}

// Dyadic green function of dipole-octupole interaction
template <typename Real>
Eigen::Matrix<Real, 3, 27> dipole_octupole_tensor(const Eigen::Vector3<Real> &unitvec) {
    Eigen::Matrix<Real, 3, 27> entries = Eigen::Matrix<Real, 3, 27>::Zero();
    for (Eigen::Index q = 0; q < 3; ++q) {
        Eigen::Index row = q;
        for (Eigen::Index j = 0; j < 3; ++j) {
            for (Eigen::Index i = 0; i < 3; ++i) {
                for (Eigen::Index k = 0; k < 3; ++k) {
                    Eigen::Index col = 9 * j + 3 * i + k;
                    entries(row, col) += dipole_octupole_coefficient(unitvec, q, j, i, k);
                }
            }
        }
    }
    return entries;
}

// Dyadic green function of octupole-dipole interaction
template <typename Real>
Eigen::Matrix<Real, 27, 3> octupole_dipole_tensor(const Eigen::Vector3<Real> &unitvec) {
    Eigen::Matrix<Real, 27, 3> entries = Eigen::Matrix<Real, 27, 3>::Zero();
    for (Eigen::Index q = 0; q < 3; ++q) {
        for (Eigen::Index j = 0; j < 3; ++j) {
            for (Eigen::Index i = 0; i < 3; ++i) {
                Eigen::Index row = 9 * q + 3 * j + i;
                for (Eigen::Index k = 0; k < 3; ++k) {
                    Eigen::Index col = k;
                    entries(row, col) += dipole_octupole_coefficient(unitvec, q, j, i, k);
                }
            }
        }
    }
    return entries;
}
} // namespace

template <typename Scalar>
GreenTensorInterpolator<Scalar> GreenTensorInterpolator<Scalar>::from_multipole_expansion(
    const std::array<real_t, 3> &distance_vector, int interaction_order) {
    GreenTensorInterpolator<Scalar> green_tensor_interpolator;

    // Normalize the distance vector, return zero green tensor if the distance is infinity
    Eigen::Map<const Eigen::Vector3<real_t>> vector_map(distance_vector.data(),
                                                        distance_vector.size());
    real_t distance = vector_map.norm();
    SPDLOG_DEBUG("Interatomic distance: {}", distance);
    if (!std::isfinite(distance)) {
        return green_tensor_interpolator;
    }
    if (distance == 0) {
        throw std::invalid_argument("The distance between the atoms must not be zero.");
    }
    Eigen::Vector3<real_t> unitvec = vector_map / distance;

    if (interaction_order >= 3) {
        green_tensor_interpolator.create_entries_from_cartesian(
            1, 1, (dipole_dipole_tensor(unitvec) / std::pow(distance, 3)).template cast<Scalar>());
    }

    if (interaction_order >= 4) {
        green_tensor_interpolator.create_entries_from_cartesian(
            1, 2,
            (dipole_quadrupole_tensor(unitvec) / std::pow(distance, 4)).template cast<Scalar>());
        green_tensor_interpolator.create_entries_from_cartesian(
            2, 1,
            (quadrupole_dipole_tensor(unitvec) / std::pow(distance, 4)).template cast<Scalar>());
    }

    if (interaction_order >= 5) {
        green_tensor_interpolator.create_entries_from_cartesian(
            2, 2,
            (quadrupole_quadrupole_tensor(unitvec) / std::pow(distance, 5))
                .template cast<Scalar>());
        green_tensor_interpolator.create_entries_from_cartesian(
            1, 3,
            (dipole_octupole_tensor(unitvec) / std::pow(distance, 5)).template cast<Scalar>());
        green_tensor_interpolator.create_entries_from_cartesian(
            3, 1,
            (octupole_dipole_tensor(unitvec) / std::pow(distance, 5)).template cast<Scalar>());
    }

    return green_tensor_interpolator;
}

// Explicit instantiations
template class GreenTensorInterpolator<double>;
template class GreenTensorInterpolator<std::complex<double>>;
} // namespace pairinteraction
