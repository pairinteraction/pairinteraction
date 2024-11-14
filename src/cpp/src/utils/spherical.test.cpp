#include "pairinteraction/utils/spherical.hpp"

#include "pairinteraction/utils/streamed.hpp"

#include <doctest/doctest.h>
#include <spdlog/spdlog.h>

namespace pairinteraction {
DOCTEST_TEST_CASE("convert cartesian to spherical basis") {
    DOCTEST_SUBCASE("kappa == 1") {
        SPDLOG_LOGGER_DEBUG(spdlog::get("doctest"), "CARTESIAN_TO_SPHERICAL_KAPPA1:\n{}",
                            fmt::streamed(spherical::CARTESIAN_TO_SPHERICAL_KAPPA1));

        auto identity = spherical::CARTESIAN_TO_SPHERICAL_KAPPA1 *
            spherical::CARTESIAN_TO_SPHERICAL_KAPPA1.adjoint();

        DOCTEST_CHECK(identity.isApprox(Eigen::Matrix3<double>::Identity(), 1e-9));
    }

    DOCTEST_SUBCASE("kappa == 2") {
        SPDLOG_LOGGER_DEBUG(spdlog::get("doctest"), "CARTESIAN_TO_SPHERICAL_KAPPA2:\n{}",
                            fmt::streamed(spherical::CARTESIAN_TO_SPHERICAL_KAPPA2));

        auto diagonal = spherical::CARTESIAN_TO_SPHERICAL_KAPPA2 *
            spherical::CARTESIAN_TO_SPHERICAL_KAPPA2.adjoint();

        DOCTEST_CHECK(diagonal.isDiagonal(1e-9));
    }
}
} // namespace pairinteraction
