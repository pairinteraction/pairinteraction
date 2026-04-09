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
namespace {

std::filesystem::path &ca_bundle_path() {
    static std::filesystem::path path;
    return path;
}

bool load_ca_bundle_from_path(httplib::SSLClient &client, const std::filesystem::path &path) {
    if (path.empty() || !std::filesystem::is_regular_file(path)) {
        return false;
    }

    std::ifstream in(path, std::ios::binary);
    std::string pem((std::istreambuf_iterator<char>(in)), {});

    if (!in || pem.empty()) {
        SPDLOG_WARN("Failed to read CA bundle from {}.", path.string());
        return false;
    }

    client.load_ca_cert_store(pem.data(), pem.size());
    SPDLOG_INFO("Using CA bundle from {}.", path.string());
    return true;
}

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

} // namespace

void set_ca_bundle_path(std::filesystem::path path) { ca_bundle_path() = std::move(path); }

std::filesystem::path get_ca_bundle_path() { return ca_bundle_path(); }

GitHubDownloader::GitHubDownloader() : client(std::make_unique<httplib::SSLClient>(host)) {
    client->set_follow_location(true);
    client->set_connection_timeout(5, 0); // seconds
    client->set_read_timeout(60, 0);      // seconds
    client->set_write_timeout(1, 0);      // seconds
    client->set_logger(log);

    if (const auto configured_path = get_ca_bundle_path(); !configured_path.empty()) {
        load_ca_bundle_from_path(*client, configured_path);
    }
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

std::string GitHubDownloader::get_host() const { return "https://" + host; }

} // namespace pairinteraction
