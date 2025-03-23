#include "pairinteraction/system/GreenTensor.hpp"

#include "pairinteraction/utils/eigen_assertion.hpp"
#include "pairinteraction/utils/traits.hpp"

#include <Eigen/Dense>
#include <complex>
#include <unsupported/Eigen/Splines>

namespace pairinteraction {
template <typename Scalar>
GreenTensor<Scalar>::Entry::Entry(int row, int col, Eigen::Spline<real_t, 1> real_spline,
                                  Eigen::Spline<real_t, 1> imag_spline)
    : row_(row), col_(col), real_spline(std::move(real_spline)),
      imag_spline(std::move(imag_spline)) {}

template <typename Scalar>
Scalar GreenTensor<Scalar>::Entry::val(double omega) const {
    if constexpr (traits::NumTraits<Scalar>::is_complex_v) {
        return {real_spline(omega)(0), imag_spline(omega)(0)};
    } else {
        return real_spline(omega)(0);
    }
}

template <typename Scalar>
int GreenTensor<Scalar>::Entry::row() const noexcept {
    return row_;
}

template <typename Scalar>
int GreenTensor<Scalar>::Entry::col() const noexcept {
    return col_;
}

template <typename Scalar>
void GreenTensor<Scalar>::set_entries(
    int kappa1, int kappa2,
    const std::vector<Eigen::SparseMatrix<Scalar>> &tensors_in_cartesian_coordinates,
    const std::vector<double> &omegas) {
    if (tensors_in_cartesian_coordinates.size() != omegas.size()) {
        throw std::invalid_argument("The number of tensors and omegas must match.");
    }

    int num_knots = static_cast<int>(omegas.size());
    Eigen::Map<const Eigen::RowVectorXd> knots(omegas.data(), num_knots);

    constexpr int spline_degree = 3; // cubic spline interpolation

    // Transform the tensors to spherical coordinates
    // TODO

    // Temporary storage wih key = (row, col) and value = vector of one double per omega
    std::map<std::pair<int, int>, std::pair<Eigen::RowVectorXd, Eigen::RowVectorXd>> temp_map;
    for (int idx = 0; idx < num_knots; ++idx) {
        const auto &tensor = tensors_in_cartesian_coordinates[idx];
        for (int k = 0; k < tensor.outerSize(); ++k) {
            for (typename Eigen::SparseMatrix<Scalar>::InnerIterator it(tensor, k); it; ++it) {
                std::pair<int, int> key{it.row(), it.col()};
                auto &vec_pair = temp_map
                                     .try_emplace(key, Eigen::RowVectorXd::Zero(num_knots),
                                                  Eigen::RowVectorXd::Zero(num_knots))
                                     .first->second;
                if constexpr (traits::NumTraits<Scalar>::is_complex_v) {
                    vec_pair.first(idx) = it.value().real();
                    vec_pair.second(idx) = it.value().imag();
                } else {
                    vec_pair.first(idx) = it.value();
                }
            }
        }
    }

    // Set the green tensor entries with spline interpolation
    std::vector<Entry> entries;
    entries.reserve(temp_map.size());
    for (const auto &[key, vec_pair] : temp_map) {
        Eigen::Spline<real_t, 1> real_spline =
            Eigen::SplineFitting<Eigen::Spline<real_t, 1>>::Interpolate(vec_pair.first,
                                                                        spline_degree, knots);

        Eigen::Spline<real_t, 1> imag_spline;
        if constexpr (traits::NumTraits<Scalar>::is_complex_v) {
            imag_spline = Eigen::SplineFitting<Eigen::Spline<real_t, 1>>::Interpolate(
                vec_pair.second, spline_degree, knots);
        }

        entries.emplace_back(key.first, key.second, std::move(real_spline), std::move(imag_spline));
    }
    entries_map[{kappa1, kappa2}] = std::move(entries);
}

template <typename Scalar>
const std::vector<typename GreenTensor<Scalar>::Entry> &
GreenTensor<Scalar>::get_entries(int kappa1, int kappa2) const {
    auto it = entries_map.find({kappa1, kappa2});
    if (it != entries_map.end()) {
        return it->second;
    }
    static const std::vector<Entry> empty_entries;
    return empty_entries;
}

// Explicit instantiations
template class GreenTensor<double>;
template class GreenTensor<std::complex<double>>;
} // namespace pairinteraction
