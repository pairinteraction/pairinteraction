// SPDX-FileCopyrightText: 2025 PairInteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include <filesystem>
#include <future>
#include <memory>
#include <string>

namespace httplib {
class SSLClient;
} // namespace httplib

namespace pairinteraction {

void set_ca_bundle_path(std::filesystem::path path);
std::filesystem::path get_ca_bundle_path();

class GitHubDownloader {
public:
    struct RateLimit {
        int remaining = -1;  // remaining number of requests
        int reset_time = -1; // unix timestamp when the rate limit resets
    };

    class Result {
    public:
        int status_code = 400;
        std::string last_modified;
        std::string body;
        RateLimit rate_limit;
    };

    GitHubDownloader();
    virtual ~GitHubDownloader();
    virtual std::future<Result> download(const std::string &remote_url,
                                         const std::string &if_modified_since = "",
                                         bool use_octet_stream = false) const;
    std::string get_host() const;

private:
    const std::string host{"api.github.com"};
    std::unique_ptr<httplib::SSLClient> client;
};

} // namespace pairinteraction
