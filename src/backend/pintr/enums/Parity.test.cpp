#include "pintr/enums/Parity.hpp"

#include "pintr/utils/streamed.hpp"

#include <algorithm>
#include <doctest/doctest.h>
#include <spdlog/spdlog.h>
#include <vector>

namespace pintr {
DOCTEST_TEST_CASE("compare parities") {
    auto p1 = Parity::ODD;
    auto p2 = Parity::EVEN;
    auto p3 = Parity::UNKNOWN;
    DOCTEST_CHECK(p1 < p2);
    DOCTEST_CHECK(p2 < p3);
}

DOCTEST_TEST_CASE("print parities") {
    auto p1 = Parity::ODD;
    auto p2 = Parity::EVEN;
    auto p3 = Parity::UNKNOWN;
    SPDLOG_LOGGER_INFO(spdlog::get("doctest"), "Parity: {}", fmt::streamed(p1));
    SPDLOG_LOGGER_INFO(spdlog::get("doctest"), "Parity: {}", fmt::streamed(p2));

    std::ostream os(nullptr);
    CHECK_THROWS_AS(os << p3, std::runtime_error);
}

DOCTEST_TEST_CASE("sort parities") {
    std::vector<Parity> parities{Parity::UNKNOWN, Parity::EVEN, Parity::ODD};
    std::sort(parities.begin(), parities.end());
    DOCTEST_CHECK(parities[0] == Parity::ODD);
    DOCTEST_CHECK(parities[1] == Parity::EVEN);
    DOCTEST_CHECK(parities[2] == Parity::UNKNOWN);
}
} // namespace pintr
