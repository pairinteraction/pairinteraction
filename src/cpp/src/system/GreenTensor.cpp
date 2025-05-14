// SPDX-FileCopyrightText: 2025 Pairinteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "pairinteraction/system/GreenTensor.hpp"

#include "pairinteraction/utils/eigen_assertion.hpp"
#include "pairinteraction/utils/eigen_compat.hpp"
#include "pairinteraction/utils/spherical.hpp"
#include "pairinteraction/utils/traits.hpp"

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <complex>
#include <map>
#include <spdlog/spdlog.h>
#include <unsupported/Eigen/Splines>

namespace pairinteraction {
template <typename Scalar>
GreenTensor<Scalar>::ConstantEntry::ConstantEntry(int row, int col, Scalar val)
    : row_(row), col_(col), val_(val) {}

template <typename Scalar>
Scalar GreenTensor<Scalar>::ConstantEntry::val() const {
    return val_;
}

template <typename Scalar>
int GreenTensor<Scalar>::ConstantEntry::row() const noexcept {
    return row_;
}

template <typename Scalar>
int GreenTensor<Scalar>::ConstantEntry::col() const noexcept {
    return col_;
}

template <typename Scalar>
GreenTensor<Scalar>::OmegaDependentEntry::OmegaDependentEntry(int row, int col,
                                                              Eigen::Spline<real_t, 1> real_spline,
                                                              Eigen::Spline<real_t, 1> imag_spline)
    : row_(row), col_(col), real_spline(std::move(real_spline)),
      imag_spline(std::move(imag_spline)) {}

template <typename Scalar>
Scalar GreenTensor<Scalar>::OmegaDependentEntry::val(double omega) const {
    if constexpr (traits::NumTraits<Scalar>::is_complex_v) {
        return {real_spline(omega)(0), imag_spline(omega)(0)};
    } else {
        return real_spline(omega)(0);
    }
}

template <typename Scalar>
int GreenTensor<Scalar>::OmegaDependentEntry::row() const noexcept {
    return row_;
}

template <typename Scalar>
int GreenTensor<Scalar>::OmegaDependentEntry::col() const noexcept {
    return col_;
}

template <typename Scalar>
void GreenTensor<Scalar>::set_entries(
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
void GreenTensor<Scalar>::set_entries(
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
const std::vector<typename GreenTensor<Scalar>::Entry> &
GreenTensor<Scalar>::get_entries(int kappa1, int kappa2) const {
    if (auto it = entries_map.find({kappa1, kappa2}); it != entries_map.end()) {
        return it->second;
    }
    static const std::vector<Entry> empty_entries;
    return empty_entries;
}

// Explicit instantiations
template class GreenTensor<double>;
template class GreenTensor<std::complex<double>>;
} // namespace pairinteraction
