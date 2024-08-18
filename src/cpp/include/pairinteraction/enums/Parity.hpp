#pragma once

#include <functional>
#include <iostream>
#include <stdexcept>

namespace pairinteraction {
enum class Parity : int { ODD = -1, EVEN = 1, UNKNOWN = 2 };

inline std::ostream &operator<<(std::ostream &os, const Parity &parity) {
    if (parity == Parity::ODD) {
        return os << "-1";
    }
    if (parity == Parity::EVEN) {
        return os << "1";
    }
    throw std::runtime_error("Cannot print unknown parity.");
}
} // namespace pairinteraction
