#ifndef KET_HPP
#define KET_HPP

#include <iostream>
#include <limits>
#include <string>

template <typename T>
class Ket {
public:
    Ket() = default;
    Ket(T energy, float f, float m, int p, std::string label);
    T get_energy();
    float get_quantum_number_f();
    float get_quantum_number_m();
    int get_parity();
    std::string get_label();

    friend std::ostream &operator<<(std::ostream &os, Ket<T> &ket) { return os << ket.get_label(); }

protected:
    virtual void ensure_assembled_ket() = 0;
    void ensure_assembled();
    void ensure_not_assembled() const;
    T energy{std::numeric_limits<T>::max()};
    float quantum_number_f{std::numeric_limits<float>::max()};
    float quantum_number_m{std::numeric_limits<float>::max()};
    int parity{std::numeric_limits<int>::max()};
    std::string label{""};

private:
    bool is_assembled{false};
};

#endif // KET_HPP
