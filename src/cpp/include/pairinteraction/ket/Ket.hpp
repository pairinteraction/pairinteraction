// SPDX-FileCopyrightText: 2024 Pairinteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include <iostream>
#include <memory>
#include <string>
#include <type_traits>

namespace pairinteraction {
enum class Parity : int;

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

    bool has_quantum_number_f() const;
    bool has_quantum_number_m() const;
    bool has_parity() const;

    double get_energy() const;
    double get_quantum_number_f() const;
    double get_quantum_number_m() const;
    Parity get_parity() const;

    virtual std::string get_label() const = 0;

    friend std::ostream &operator<<(std::ostream &os, const Ket &ket) {
        return os << ket.get_label();
    }

protected:
    Ket(double energy, double f, double m, Parity p);

    bool operator==(const Ket &other) const;

    struct hash {
        std::size_t operator()(const Ket &k) const;
    };

    double energy;
    double quantum_number_f;
    double quantum_number_m;
    Parity parity;
};
} // namespace pairinteraction
