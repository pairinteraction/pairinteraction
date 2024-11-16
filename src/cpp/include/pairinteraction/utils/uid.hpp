#pragma once

#include <chrono>
#include <cstdint>
#include <fmt/core.h>
#include <random>
#include <string>

namespace pairinteraction {
namespace utils {

inline std::string generate_uid() {
    // Get a 64-bit timestamp
    auto now = std::chrono::high_resolution_clock::now();
    std::uint64_t timestamp =
        std::chrono::duration_cast<std::chrono::nanoseconds>(now.time_since_epoch()).count();

    // Get a 64-bit random number
    static std::mt19937_64 rng(std::random_device{}());
    std::uniform_int_distribution<uint64_t> dis(0, UINT64_MAX);
    std::uint64_t random = dis(rng);

    return fmt::format("{:016x}{:016x}", timestamp, random);
}

} // namespace utils
} // namespace pairinteraction
