// SPDX-FileCopyrightText: 2026 PairInteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include <stdexcept>
#include <string>
#include <utility>

namespace pairinteraction {

/**
 * @class QuantumNumberNotAvailableError
 *
 * @brief Exception thrown when a species does not provide a requested quantum number.
 *
 * The exception carries the name of the requested quantum number and the species so that the
 * caller (e.g. the Python bindings) can react to the specific quantum number that is missing.
 */
class QuantumNumberNotAvailableError : public std::invalid_argument {
public:
    QuantumNumberNotAvailableError(std::string name, std::string species)
        : std::invalid_argument("The quantum number '" + name +
                                "' is not available for the species '" + species +
                                "', because the table of states of the species does not provide "
                                "it."),
          name(std::move(name)), species(std::move(species)) {}

    const std::string &get_name() const { return name; }
    const std::string &get_species() const { return species; }

private:
    std::string name;
    std::string species;
};

} // namespace pairinteraction
