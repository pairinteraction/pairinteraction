// SPDX-FileCopyrightText: 2024 PairInteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include <memory>
#include <string>
#include <type_traits>
#include <unordered_map>

namespace pairinteraction {

/**
 * @class Ket
 *
 * @brief Base class for a ket.
 *
 * This base class represents a ket. It is a base class for specific ket implementations. Its
 * constructor is protected to indicate that derived classes should not allow direct instantiation.
 * Instead, a factory class should be provided that is a friend of the derived class and can create
 * instances of it.
 *
 * The ket stores the quantum numbers that are available for it in a map.
 */

class Ket {
public:
    Ket() = delete;
    virtual ~Ket() = default;

    double get_energy() const;
    bool has_quantum_number(const std::string &name) const;
    double get_quantum_number(const std::string &name) const;

protected:
    Ket(double energy, std::unordered_map<std::string, double> quantum_numbers);

    bool operator==(const Ket &other) const;

    struct hash {
        std::size_t operator()(const Ket &k) const;
    };

    double energy;
    std::unordered_map<std::string, double> quantum_numbers;
};
} // namespace pairinteraction
