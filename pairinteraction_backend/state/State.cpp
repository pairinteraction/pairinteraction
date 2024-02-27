#include "State.hpp"

template <typename T>
State<T>::State(T energy, float f, float m, int p, std::string label)
    : energy(energy), quantum_number_f(f), quantum_number_m(m), parity(p), label(label) {}

template <typename T>
T State<T>::get_energy() const {
    return energy;
}

template <typename T>
float State<T>::get_quantum_number_f() const {
    return quantum_number_f;
}

template <typename T>
float State<T>::get_quantum_number_m() const {
    return quantum_number_m;
}

template <typename T>
int State<T>::get_parity() const {
    return parity;
}

template <typename T>
std::string State<T>::get_label() const {
    return label;
}

// Explicit instantiations
template class State<float>;
template class State<double>;

///////////////////////////////////////////////////////////////////////////////////////
// Test cases
///////////////////////////////////////////////////////////////////////////////////////

#include <doctest/doctest.h>

DOCTEST_TEST_CASE("constructing a class derived from state") {
    class StateDerived : public State<float> {
    public:
        StateDerived(float energy, float f, float m, int p, std::string label)
            : State<float>(energy, f, m, p, label) {}
    };

    StateDerived state(1.0f, 2.0f, 3.0f, 4, "label");
}
