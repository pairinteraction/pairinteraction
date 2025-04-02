// SPDX-FileCopyrightText: 2024 Pairinteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include <string>
#include <type_traits>

namespace pairinteraction::utils {
constexpr size_t OFFSET = 500;

// String that can be used in an SQL query to calculate a linearized index ("id_in_database")
// from the id of a state without quantum number m and the quantum number m. Note that the
// string must be of the form "id*{2*OFFSET}+(2*m+{OFFSET})::bigint".
constexpr std::string_view SQL_TERM_FOR_LINEARIZED_ID_IN_DATABASE = "id*1000+(2*m+500)::bigint";

// Method for calculating the linearized index
inline size_t get_linearized_id_in_database(size_t id, double quantum_number_m) {
    return id * 2 * OFFSET + static_cast<size_t>(2 * quantum_number_m + OFFSET);
}
} // namespace pairinteraction::utils
