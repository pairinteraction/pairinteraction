#pragma once

#include "utils/traits.hpp"

#include <array>
#include <limits>
#include <stdexcept>

namespace spherical {

template <typename Scalar>
inline std::array<Scalar, 3>
convert_to_spherical_basis(const std::array<typename traits::NumTraits<Scalar>::real_t, 3> &field) {
    using real_t = typename traits::NumTraits<Scalar>::real_t;

    if constexpr (traits::NumTraits<Scalar>::is_complex_v) {
        return {field[2],
                std::complex<real_t>(-field[0], -field[1]) / static_cast<real_t>(std::sqrt(2)),
                std::complex<real_t>(field[0], -field[1]) / static_cast<real_t>(std::sqrt(2))};
    } else {
        if (std::abs(field[1]) > 10 * std::numeric_limits<real_t>::epsilon()) {
            throw std::invalid_argument(
                "The field must not have a y-component if the scalar type is real.");
        }
        return {field[2], -field[0] / static_cast<real_t>(std::sqrt(2)),
                field[0] / static_cast<real_t>(std::sqrt(2))};
    }
}

} // namespace spherical
