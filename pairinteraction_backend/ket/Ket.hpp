#pragma once

#include <iostream>
#include <string>

template <typename T>
class Ket {
public:
    Ket() = delete;
    virtual ~Ket() = default;
    T get_energy() const;
    float get_quantum_number_f() const;
    float get_quantum_number_m() const;
    int get_parity() const;
    std::string get_label() const;

    friend std::ostream &operator<<(std::ostream &os, const Ket<T> &ket) {
        return os << ket.get_label();
    }

protected:
    Ket(T energy, float f, float m, int p, std::string label);
    T energy;
    float quantum_number_f;
    float quantum_number_m;
    int parity;
    std::string label;
};

extern template class Ket<float>;
extern template class Ket<double>;
