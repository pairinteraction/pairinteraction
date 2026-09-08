// SPDX-FileCopyrightText: 2024 PairInteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "pairinteraction/utils/spherical.hpp"

#include <doctest/doctest.h>

namespace pairinteraction {
DOCTEST_TEST_CASE("convert cartesian to spherical basis") {
    DOCTEST_SUBCASE("kappa == 0") {
        // The monopole is a scalar, the transformation is trivial
        const auto &real_mat = spherical::get_transformator<double>(0);
        const auto &complex_mat = spherical::get_transformator<std::complex<double>>(0);

        DOCTEST_CHECK(real_mat.rows() == 1);
        DOCTEST_CHECK(real_mat.cols() == 1);
        DOCTEST_CHECK(real_mat(0, 0) == 1);
        DOCTEST_CHECK(complex_mat.rows() == 1);
        DOCTEST_CHECK(complex_mat.cols() == 1);
        DOCTEST_CHECK(complex_mat(0, 0) == std::complex<double>(1, 0));
    }

    DOCTEST_SUBCASE("kappa == 1") {
        auto identity = spherical::CARTESIAN_TO_SPHERICAL_KAPPA1 *
            spherical::CARTESIAN_TO_SPHERICAL_KAPPA1.adjoint();

        DOCTEST_CHECK(identity.isApprox(Eigen::Matrix3<double>::Identity(), 1e-9));
    }

    DOCTEST_SUBCASE("kappa == 2") {
        auto diagonal = spherical::CARTESIAN_TO_SPHERICAL_KAPPA2 *
            spherical::CARTESIAN_TO_SPHERICAL_KAPPA2.adjoint();

        DOCTEST_CHECK(diagonal.isDiagonal(1e-9));
    }
}
} // namespace pairinteraction
