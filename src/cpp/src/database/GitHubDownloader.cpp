// SPDX-FileCopyrightText: 2025 PairInteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "pairinteraction/database/GitHubDownloader.hpp"

#include "pairinteraction/utils/paths.hpp"

#include <filesystem>
#include <fmt/core.h>
#include <fstream>
#include <future>
#include <httplib.h>
#include <spdlog/spdlog.h>
#include <stdexcept>

namespace pairinteraction {

void log(const httplib::Request &req, const httplib::Response &res) {
    if (!spdlog::default_logger()->should_log(spdlog::level::debug)) {
        return;
    }

    SPDLOG_DEBUG("[httplib] {} {}", req.method, req.path);
    for (const auto &[k, v] : req.headers) {
        SPDLOG_DEBUG("[httplib]   {}: {}\n", k, v);
    }

    SPDLOG_DEBUG("[httplib] Response with status {}", res.status);
    for (const auto &[k, v] : res.headers) {
        SPDLOG_DEBUG("[httplib]   {}: {}\n", k, v);
    }

    if (res.body.empty()) {
        return;
    }

    if (res.body.size() > 1024) {
        SPDLOG_DEBUG("[httplib] Body ({} bytes, first 1024 bytes):", res.body.size());
    } else {
        SPDLOG_DEBUG("[httplib] Body ({} bytes):", res.body.size());
    }
    SPDLOG_DEBUG("[httplib]   {}", res.body.substr(0, 1024));
}

GitHubDownloader::GitHubDownloader() : client(std::make_unique<httplib::SSLClient>(host)) {
    std::filesystem::path configdir = paths::get_config_directory();
    if (!std::filesystem::exists(configdir)) {
        std::filesystem::create_directories(configdir);
    } else if (!std::filesystem::is_directory(configdir)) {
        throw std::filesystem::filesystem_error("Cannot access config directory ",
                                                configdir.string(),
                                                std::make_error_code(std::errc::not_a_directory));
    }

    std::filesystem::path cert_path = configdir / "ca-bundle.crt";
    if (!std::filesystem::exists(cert_path)) {
        std::ofstream out(cert_path);
        if (!out) {
            throw std::runtime_error("Failed to create certificate file at " + cert_path.string());
        }
        out << cert;
        out.close();
    }
    cert_path_str = cert_path.string();

    client->set_follow_location(true);
    client->set_connection_timeout(5, 0); // seconds
    client->set_read_timeout(60, 0);      // seconds
    client->set_write_timeout(1, 0);      // seconds
    client->set_ca_cert_path(cert_path_str);
    client->set_logger(log);
}

GitHubDownloader::~GitHubDownloader() = default;

std::future<GitHubDownloader::Result>
GitHubDownloader::download(const std::string &remote_url, const std::string &if_modified_since,
                           bool use_octet_stream) const {
    return std::async(
        std::launch::async, [this, remote_url, if_modified_since, use_octet_stream]() -> Result {
            SPDLOG_DEBUG("Downloading from GitHub: {}", remote_url);

            // Prepare headers
            httplib::Headers headers{
                {"User-Agent", "pairinteraction"},
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

            // If we're fetching binary, stream with a progress callback; otherwise use a simple get
            httplib::Result response;
            std::string streamed_body;
            if (use_octet_stream) {
                auto content_receiver = [&](const char *data, size_t len) {
                    streamed_body.append(data, len);
                    return true;
                };

                // Progress display
                int last_pct = -1;
                auto progress_display = [&last_pct, remote_url](uint64_t cur, uint64_t total) {
                    if (total == 0) {
                        fmt::print(stderr, "\rDownloading {}...", remote_url);
                        (void)std::fflush(stderr);
                    } else if (int pct = static_cast<int>((cur * 100) / total); pct != last_pct) {
                        last_pct = pct;
                        fmt::print(stderr, "\rDownloading {}... {:3d}%", remote_url, pct);
                        (void)std::fflush(stderr);
                    }
                    return true;
                };

                response = client->Get(remote_url, headers, content_receiver, progress_display);

                // Ensure the progress display ends cleanly if we showed it
                if (last_pct >= 0) {
                    fmt::print(stderr, "\n");
                    (void)std::fflush(stderr);
                }
            } else {
                response = client->Get(remote_url, headers);
            }

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
            result.body = use_octet_stream ? std::move(streamed_body) : response->body;
            result.status_code = response->status;

            SPDLOG_DEBUG("Response status: {}", response->status);
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

std::string GitHubDownloader::get_host() const { return "https://" + host; }

} // namespace pairinteraction
