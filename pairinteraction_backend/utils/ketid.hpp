#pragma once

#include <string>
#include <type_traits>

namespace ketid::atom {
constexpr size_t OFFSET = 500;

// String that can be used in an SQL query to calculate a linearized index ("ketid") from
// the id of a state without quantum number m and the quantum number m. Note that the
// string must be of the form "id*{2*OFFSET}+(2*m+{OFFSET})::bigint".
constexpr std::string_view SQL_TERM = "id*1000+(2*m+500)::bigint";

// Method for calculating the linearized index
template <typename Real>
inline size_t get(size_t id, Real quantum_number_m) {
    static_assert(std::is_arithmetic_v<Real>);

    return id * 2 * OFFSET + static_cast<size_t>(2 * quantum_number_m + OFFSET);
}

// Method to update the linearized index when the quantum number m changes
template <typename Real>
inline size_t transform(size_t ketid, Real old_quantum_number_m, Real new_quantum_number_m) {
    static_assert(std::is_arithmetic_v<Real>);

    return ketid + static_cast<size_t>(2 * new_quantum_number_m + OFFSET) -
        static_cast<size_t>(2 * old_quantum_number_m + OFFSET);
}
} // namespace ketid::atom
