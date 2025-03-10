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
    RateLimit get_rate_limit() const;
    std::string get_host() const;

private:
    const std::string host{"https://api.github.com"};
    const std::string github_ca_cert{R"(-----BEGIN CERTIFICATE-----
MIIEoDCCBEagAwIBAgIQKhb1wgEYB/cKkmPdPDmp8jAKBggqhkjOPQQDAjCBjzEL
MAkGA1UEBhMCR0IxGzAZBgNVBAgTEkdyZWF0ZXIgTWFuY2hlc3RlcjEQMA4GA1UE
BxMHU2FsZm9yZDEYMBYGA1UEChMPU2VjdGlnbyBMaW1pdGVkMTcwNQYDVQQDEy5T
ZWN0aWdvIEVDQyBEb21haW4gVmFsaWRhdGlvbiBTZWN1cmUgU2VydmVyIENBMB4X
DTI1MDIwNTAwMDAwMFoXDTI2MDIwNTIzNTk1OVowFzEVMBMGA1UEAwwMKi5naXRo
dWIuY29tMFkwEwYHKoZIzj0CAQYIKoZIzj0DAQcDQgAElbQ+DErBU9/BYlYXV5qx
aS5Nu5Ucd+scwICjYp7Z2YcJ2Jgu7HBjM1R6+b2d8yZYJlpLB2aX3qGwc1ZscE7w
HKOCAvkwggL1MB8GA1UdIwQYMBaAFPaFCjsRhuEEfQ6qCyzS7sxke3uuMB0GA1Ud
DgQWBBSY8zci2pdh0RUe+FTOokexsSdvADAOBgNVHQ8BAf8EBAMCB4AwDAYDVR0T
AQH/BAIwADAdBgNVHSUEFjAUBggrBgEFBQcDAQYIKwYBBQUHAwIwSQYDVR0gBEIw
QDA0BgsrBgEEAbIxAQICBzAlMCMGCCsGAQUFBwIBFhdodHRwczovL3NlY3RpZ28u
Y29tL0NQUzAIBgZngQwBAgEwgYQGCCsGAQUFBwEBBHgwdjBPBggrBgEFBQcwAoZD
aHR0cDovL2NydC5zZWN0aWdvLmNvbS9TZWN0aWdvRUNDRG9tYWluVmFsaWRhdGlv
blNlY3VyZVNlcnZlckNBLmNydDAjBggrBgEFBQcwAYYXaHR0cDovL29jc3Auc2Vj
dGlnby5jb20wggF9BgorBgEEAdZ5AgQCBIIBbQSCAWkBZwB1AJaXZL9VWJet90OH
aDcIQnfp8DrV9qTzNm5GpD8PyqnGAAABlNNtegoAAAQDAEYwRAIgARjCCRJjeSsj
1aDf1k2e9k7BMQOzBk5NH2yNZnExXJYCIEx3XBFFZ4kmjLneOJG0C/W/K+1W+Lx5
vY+UM1oFAV7SAHYAGYbUxyiqb/66A294Kk0BkarOLXIxD67OXXBBLSVMx9QAAAGU
0215mgAABAMARzBFAiA9Ui4Z4WE8Mg11ZmjscGtczwGRFgRySEgdY3O/Hmn3rgIh
AJSNm0C0lHX7yf3IW+dMIctrApWm2rL2P6Wei1wEBuovAHYAyzj3FYl8hKFEX1vB
3fvJbvKaWc1HCmkFhbDLFMMUWOcAAAGU0215zwAABAMARzBFAiBNfjmh0gVRA7aq
sfUsyIvg9rvkxJjrV0/w4sjsTbVqngIhAM8fdljMks2Hh+It5rYLMS4utsRzTDEU
8zBOHS3K9bSDMCMGA1UdEQQcMBqCDCouZ2l0aHViLmNvbYIKZ2l0aHViLmNvbTAK
BggqhkjOPQQDAgNIADBFAiAHU6XJYeE/cjpdM9Rfn0IdGZEEs0zTAuyN0mUkXnZa
KQIhAPkzw1jn3t/HCIVT98ZYrTYNsxElJRY6JH89pJlVt3Ay
-----END CERTIFICATE-----)"};
    std::unique_ptr<httplib::Client> client;
};

} // namespace pairinteraction
