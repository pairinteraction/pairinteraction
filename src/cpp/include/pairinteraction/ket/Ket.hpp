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
 *
 * @tparam Real Real number type.
 */

template <typename Real>
class Ket {
    static_assert(std::is_floating_point_v<Real>);

public:
    Ket() = delete;
    virtual ~Ket() = default;

    bool has_quantum_number_f() const;
    bool has_quantum_number_m() const;
    bool has_parity() const;

    Real get_energy() const;
    Real get_quantum_number_f() const;
    Real get_quantum_number_m() const;
    Parity get_parity() const;

    virtual std::string get_label() const = 0;

    friend std::ostream &operator<<(std::ostream &os, const Ket<Real> &ket) {
        return os << ket.get_label();
    }

protected:
    Ket(Real energy, Real f, Real m, Parity p);

    bool operator==(const Ket<Real> &other) const;

    struct hash {
        std::size_t operator()(const Ket<Real> &k) const;
    };

    Real energy;
    Real quantum_number_f;
    Real quantum_number_m;
    Parity parity;
};

extern template class Ket<double>;
} // namespace pairinteraction
