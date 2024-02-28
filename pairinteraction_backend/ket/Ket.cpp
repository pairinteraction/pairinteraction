#include "Ket.hpp"

#include <limits>

template <typename T>
Ket<T>::Ket(T energy, float f, float m, int p, std::string label)
    : energy(energy), quantum_number_f(f), quantum_number_m(m), parity(p), label(label),
      is_assembled(true) {}

template <typename T>
void Ket<T>::ensure_assembled() {
    if (!is_assembled) {
        this->ensure_assembled_ket();
    }
}

template <typename T>
void Ket<T>::ensure_not_assembled() const {
    if (is_assembled) {
        throw std::runtime_error("Ket is already assembled");
    }
}

template <typename T>
T Ket<T>::get_energy() {
    this->ensure_assembled();
    return energy;
}

template <typename T>
float Ket<T>::get_quantum_number_f() {
    this->ensure_assembled();
    return quantum_number_f;
}

template <typename T>
float Ket<T>::get_quantum_number_m() {
    this->ensure_assembled();
    return quantum_number_m;
}

template <typename T>
int Ket<T>::get_parity() {
    this->ensure_assembled();
    return parity;
}

template <typename T>
std::string Ket<T>::get_label() {
    this->ensure_assembled();
    return label;
}

// Explicit instantiations
template class Ket<float>;
template class Ket<double>;

///////////////////////////////////////////////////////////////////////////////////////
// Test cases
///////////////////////////////////////////////////////////////////////////////////////

#include <doctest/doctest.h>
#include <sstream>

DOCTEST_TEST_CASE("constructing a class derived from ket") {
    class KetDerived : public Ket<float> {
    public:
        KetDerived(float energy, float f, float m, int p, std::string label)
            : Ket<float>(energy, f, m, p, label) {}

    private:
        void ensure_assembled_ket() override {}
    };

    KetDerived ket(1.0f, 2.0f, 3.0f, 4, "my_label");

    // Check that the label can be printed
    std::stringstream ss;
    ss << ket;
    CHECK(ss.str() == "my_label");
}
