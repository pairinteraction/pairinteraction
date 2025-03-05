#include "pairinteraction/database/GitHubDownloader.hpp"

#include <fmt/core.h>
#include <future>
#include <httplib.h>
#include <stdexcept>

namespace pairinteraction {

GitHubDownloader::GitHubDownloader() : client(std::make_unique<httplib::Client>(host)) {
    client->set_follow_location(true);
    client->set_connection_timeout(5, 0); // seconds
    client->set_read_timeout(60, 0);      // seconds
    client->set_write_timeout(1, 0);      // seconds
    client->enable_server_certificate_verification(true);
}

GitHubDownloader::~GitHubDownloader() = default;

std::future<GitHubDownloader::Result>
GitHubDownloader::download(const std::string &remote_url, const std::string &if_modified_since,
                           bool use_octet_stream) const {
    return std::async(
        std::launch::async, [this, remote_url, if_modified_since, use_octet_stream]() -> Result {
            // Prepare headers
            httplib::Headers headers{
                {"X-GitHub-Api-Version", "2022-11-28"},
                {"Accept",
                 use_octet_stream ? "application/octet-stream" : "application/vnd.github+json"}};

            if (!if_modified_since.empty()) {
                headers.emplace("if-modified-since", if_modified_since);
            }

            // Use the GitHub token if available; otherwise, if we have a conditional request,
            // insert a dummy authorization header to avoid increasing rate limits
            if (auto *token = std::getenv("GITHUB_TOKEN"); token) {
                headers.emplace("Authorization", fmt::format("Bearer {}", token));
            } else if (!if_modified_since.empty()) {
                headers.emplace("Authorization",
                                "avoids-an-increase-in-ratelimits-used-if-304-is-returned");
            }

            auto response = client->Get(remote_url.c_str(), headers);

            // Handle if the response is null
            if (!response) {
                // Defensive handling: if response is null and the error is unknown,
                // treat this as a 304 Not Modified
                if (response.error() == httplib::Error::Unknown) {
                    return Result{304, "", "", {}};
                }
                throw std::runtime_error(fmt::format("Error downloading '{}': {}", remote_url,
                                                     httplib::to_string(response.error())));
            }

            // Parse the response
            Result result;
            if (response->has_header("x-ratelimit-remaining")) {
                result.rate_limit.remaining =
                    std::stoi(response->get_header_value("x-ratelimit-remaining"));
            }
            if (response->has_header("x-ratelimit-reset")) {
                result.rate_limit.reset_time =
                    std::stoi(response->get_header_value("x-ratelimit-reset"));
            }
            if (response->has_header("last-modified")) {
                result.last_modified = response->get_header_value("last-modified");
            }
            result.body = response->body;
            result.status_code = response->status;
            return result;
        });
}

GitHubDownloader::RateLimit GitHubDownloader::get_rate_limit() const {
    // This call now either returns valid rate limit data or throws an exception on error
    Result result = download("/rate_limit", "", false).get();
    if (result.status_code != 200) {
        throw std::runtime_error(
            fmt::format("Failed obtaining the rate limit: status code {}.", result.status_code));
    }
    return result.rate_limit;
}

std::string GitHubDownloader::get_host() const { return host; }

} // namespace pairinteraction
