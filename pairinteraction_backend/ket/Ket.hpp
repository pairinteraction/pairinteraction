#pragma once

#include <iostream>
#include <string>

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
public:
    Ket() = delete;
    virtual ~Ket() = default;
    Real get_energy() const;
    float get_quantum_number_f() const;
    float get_quantum_number_m() const;
    int get_parity() const;
    std::string get_label() const;

    friend std::ostream &operator<<(std::ostream &os, const Ket<Real> &ket) {
        return os << ket.get_label();
    }

protected:
    Ket(Real energy, float f, float m, int p, std::string label);
    Real energy;
    float quantum_number_f;
    float quantum_number_m;
    int parity;
    std::string label;
};

extern template class Ket<float>;
extern template class Ket<double>;
