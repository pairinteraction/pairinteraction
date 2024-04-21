#pragma once

#include <fmt/format.h>

namespace ketid {
namespace atom {
const size_t OFFSET = 500;
const std::string SQL_TERM = fmt::format("id*{}+(2*m+{})::bigint", 2 * OFFSET, OFFSET);

inline size_t get(size_t id, float quantum_number_m) {
    return id * 2 * OFFSET + static_cast<size_t>(2 * quantum_number_m + OFFSET);
}

inline size_t transform(size_t ketid, float old_quantum_number_m, float new_quantum_number_m) {
    return ketid + static_cast<size_t>(2 * new_quantum_number_m + OFFSET) -
        static_cast<size_t>(2 * old_quantum_number_m + OFFSET);
}
} // namespace atom
} // namespace ketid
