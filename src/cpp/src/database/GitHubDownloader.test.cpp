#include "pairinteraction/database/GitHubDownloader.hpp"

#include <ctime>
#include <doctest/doctest.h>
#include <httplib.h>

namespace pairinteraction {
DOCTEST_TEST_CASE(
    "GitHubDownloader basic functionality" *
    doctest::skip(true)) { // TODO skip if !Database::get_global_instance().download_missing()
    GitHubDownloader downloader;

    int current_time = static_cast<int>(std::time(nullptr));

    DOCTEST_SUBCASE("Successful download") {
        auto future = downloader.download("/octocat");
        auto result = future.get();

        DOCTEST_CHECK(result.success == true);
        DOCTEST_CHECK(result.status_code == 200);
        DOCTEST_CHECK(result.rate_limit.reset_time > current_time);
        DOCTEST_CHECK(result.rate_limit.remaining >= 0);
        DOCTEST_MESSAGE("Number of remaining requests: ", result.rate_limit.remaining);
    }

    DOCTEST_SUBCASE("Rate limit") {
        auto rate_limit = downloader.get_rate_limit();

        DOCTEST_CHECK(rate_limit.reset_time > current_time);
        DOCTEST_CHECK(rate_limit.remaining >= 0);
        DOCTEST_MESSAGE("Number of remaining requests: ", rate_limit.remaining);
    }
}
} // namespace pairinteraction
