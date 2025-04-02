// SPDX-FileCopyrightText: 2024 Pairinteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "pairinteraction/utils/spherical.hpp"

#include <doctest/doctest.h>

namespace pairinteraction {
DOCTEST_TEST_CASE("convert cartesian to spherical basis") {
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
