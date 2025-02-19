#pragma once

#include <future>
#include <memory>
#include <string>

namespace httplib {
class Client;
}

namespace pairinteraction {

class GitHubDownloader {
public:
    struct RateLimit {
        int remaining = -1;  // remaining number of requests
        int reset_time = -1; // unix timestamp when the rate limit resets
    };

    class Result {
    public:
        bool success;
        int status_code;
        std::string last_modified;
        std::string body;
        RateLimit rate_limit;
    };

public:
    GitHubDownloader();
    virtual ~GitHubDownloader();
    virtual std::future<Result> download(const std::string &remote_url,
                                         const std::string &if_modified_since = "",
                                         bool use_octet_stream = false);
    RateLimit get_rate_limit();
    std::string get_host() const;

private:
    std::string host{"https://api.github.com"};
    std::unique_ptr<httplib::Client> client;
};

} // namespace pairinteraction
