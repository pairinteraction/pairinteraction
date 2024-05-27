#pragma once

#include <fmt/format.h>

namespace ketid {
namespace atom {
const size_t OFFSET = 500;

// String that can be used in an SQL query to calculate a linearized index ("ketid") from
// the id of a state without quantum number m and the quantum number m
const std::string SQL_TERM = fmt::format("id*{}+(2*m+{})::bigint", 2 * OFFSET, OFFSET);

// Method for calculating the linearized index
template <typename Real>
inline size_t get(size_t id, Real quantum_number_m) {
    return id * 2 * OFFSET + static_cast<size_t>(2 * quantum_number_m + OFFSET);
}

// Method to update the linearized index when the quantum number m changes
template <typename Real>
inline size_t transform(size_t ketid, Real old_quantum_number_m, Real new_quantum_number_m) {
    return ketid + static_cast<size_t>(2 * new_quantum_number_m + OFFSET) -
        static_cast<size_t>(2 * old_quantum_number_m + OFFSET);
}
} // namespace atom
} // namespace ketid
