// SPDX-FileCopyrightText: 2024 PairInteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include <iostream>
#include <memory>
#include <string>
#include <type_traits>

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
 */

class Ket {
public:
    Ket() = delete;
    virtual ~Ket() = default;

    double get_energy() const;

    virtual std::string get_label() const = 0;

    friend std::ostream &operator<<(std::ostream &os, const Ket &ket) {
        return os << ket.get_label();
    }

protected:
    explicit Ket(double energy);

    bool operator==(const Ket &other) const;

    struct hash {
        std::size_t operator()(const Ket &k) const;
    };

    double energy;
};
} // namespace pairinteraction
