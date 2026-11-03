// SPDX-FileCopyrightText: 2024 PairInteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include <array>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>

namespace pairinteraction {
// The quantum numbers by which the states of a basis can be labeled, sorted and grouped into
// blocks. To add a state label, add an enumerator and its name to state_label_names.
enum class StateLabel : unsigned char {
    QUANTUM_NUMBER_F,
    QUANTUM_NUMBER_M,
    PARITY,
    PRODUCT_OF_PARITIES,
    PARITY_UNDER_INVERSION,
    PARITY_UNDER_PERMUTATION
};

// The names of the state labels, by which kets and creators refer to the quantum numbers
inline constexpr auto state_label_names = std::to_array<std::pair<StateLabel, std::string_view>>({
    {StateLabel::QUANTUM_NUMBER_F, "f"},
    {StateLabel::QUANTUM_NUMBER_M, "m"},
    {StateLabel::PARITY, "parity"},
    {StateLabel::PRODUCT_OF_PARITIES, "product_of_parities"},
    {StateLabel::PARITY_UNDER_INVERSION, "parity_under_inversion"},
    {StateLabel::PARITY_UNDER_PERMUTATION, "parity_under_permutation"},
});

constexpr std::string_view get_state_label_name(StateLabel label) {
    for (const auto &[state_label, name] : state_label_names) {
        if (state_label == label) {
            return name;
        }
    }
    throw std::invalid_argument("The state label " + std::to_string(std::to_underlying(label)) +
                                " has no name.");
}

// The label of the quantity with the given name
constexpr StateLabel get_state_label(std::string_view name) {
    for (const auto &[state_label, state_label_name] : state_label_names) {
        if (state_label_name == name) {
            return state_label;
        }
    }
    throw std::invalid_argument("The states cannot be labeled by " + std::string(name) + ".");
}
} // namespace pairinteraction
