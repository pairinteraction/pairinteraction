#include "pairinteraction/utils/spherical.hpp"

#include "pairinteraction/utils/streamed.hpp"

#include <doctest/doctest.h>
#include <spdlog/spdlog.h>

namespace pairinteraction {
DOCTEST_TEST_CASE("convert cartesian to spherical basis") {
    DOCTEST_SUBCASE("kappa == 1") {
        SPDLOG_LOGGER_DEBUG(spdlog::get("doctest"), "CARTESIAN_TO_SPHERICAL1:\n{}",
                            fmt::streamed(spherical::CARTESIAN_TO_SPHERICAL1));

        auto identity =
            spherical::CARTESIAN_TO_SPHERICAL1 * spherical::CARTESIAN_TO_SPHERICAL1.adjoint();
        DOCTEST_CHECK(identity.isApprox(Eigen::Matrix3<double>::Identity(), 1e-9));
    }

    DOCTEST_SUBCASE("kappa == 2") {
        SPDLOG_LOGGER_DEBUG(spdlog::get("doctest"), "CARTESIAN_TO_SPHERICAL2:\n{}",
                            fmt::streamed(spherical::CARTESIAN_TO_SPHERICAL2));

        auto identity = std::sqrt(2. / 3.) * spherical::CARTESIAN_TO_SPHERICAL2 *
            std::sqrt(2. / 3.) * spherical::CARTESIAN_TO_SPHERICAL2.adjoint();
        DOCTEST_CHECK(identity.isApprox(Eigen::Matrix<double, 5, 5>::Identity(), 1e-9));
    }
}
} // namespace pairinteraction
