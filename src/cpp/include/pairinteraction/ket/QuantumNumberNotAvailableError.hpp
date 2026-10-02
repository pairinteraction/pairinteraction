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
 * @brief Exception thrown when a ket does not provide a requested quantum number.
 *
 * The exception carries the name of the requested quantum number so that the caller (e.g. the
 * Python bindings) can react to the specific quantum number that is missing.
 */
class QuantumNumberNotAvailableError : public std::invalid_argument {
public:
    explicit QuantumNumberNotAvailableError(std::string name)
        : std::invalid_argument("The quantum number '" + name + "' is not available for the ket."),
          name(std::move(name)) {}

    const std::string &get_name() const { return name; }

private:
    std::string name;
};

} // namespace pairinteraction
