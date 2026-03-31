// SPDX-FileCopyrightText: 2025 PairInteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <array>
#include <cstdio>
#include <filesystem>
#include <httplib.h>
#include <stdexcept>
#include <string>

int main(int /*argc*/, char ** /*argv*/) {
    const std::filesystem::path cert_path{CERTIFI_CA_BUNDLE_PATH};

    // Use the certificate to make a request
    httplib::SSLClient cli("api.github.com");
    cli.set_ca_cert_path(cert_path.string());
    const auto res = cli.Get("/rate_limit");

    // Check the response
    if (!res || res->status != 200) {
        throw std::runtime_error("Error: " + httplib::to_string(res.error()));
    }
    return 0;
}
