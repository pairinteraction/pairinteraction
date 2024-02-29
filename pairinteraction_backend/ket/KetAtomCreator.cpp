#include "ket/KetAtomCreator.hpp"

#include <cmath>
#include <limits>

template <typename T>
KetAtomCreator<T>::KetAtomCreator(std::string species) : species(species) {}

template <typename T>
KetAtomCreator<T>::KetAtomCreator(std::string species, int n, T l, float j, float m)
    : species(species), quantum_number_n(n), quantum_number_l(l), quantum_number_f(j),
      quantum_number_j(j), quantum_number_m(m), quantum_number_s(0.5) {}

template <typename T>
void KetAtomCreator<T>::set_energy(T value) {
    energy.emplace(value);
}

template <typename T>
void KetAtomCreator<T>::set_quantum_number_f(float value) {
    quantum_number_f.emplace(value);
}

template <typename T>
void KetAtomCreator<T>::set_quantum_number_m(float value) {
    quantum_number_m.emplace(value);
}

template <typename T>
void KetAtomCreator<T>::set_parity(int value) {
    parity.emplace(value);
}

template <typename T>
void KetAtomCreator<T>::set_quantum_number_n(int value) {
    quantum_number_n.emplace(value);
}

template <typename T>
void KetAtomCreator<T>::set_quantum_number_nu(T value) {
    quantum_number_nu.emplace(value);
}

template <typename T>
void KetAtomCreator<T>::set_quantum_number_l(T value) {
    quantum_number_l.emplace(value);
}

template <typename T>
void KetAtomCreator<T>::set_quantum_number_s(T value) {
    quantum_number_s.emplace(value);
}

template <typename T>
void KetAtomCreator<T>::set_quantum_number_j(T value) {
    quantum_number_j.emplace(value);
}

template <typename T>
KetAtom<T> KetAtomCreator<T>::create() const {

    // TODO perform database request

    return KetAtom<T>(energy.value_or(std::numeric_limits<T>::quiet_NaN()),
                      quantum_number_f.value_or(std::numeric_limits<float>::quiet_NaN()),
                      quantum_number_m.value_or(std::numeric_limits<float>::quiet_NaN()),
                      parity.value_or(std::pow(
                          -1, quantum_number_l.value_or(std::numeric_limits<T>::quiet_NaN()))),
                      "", species, quantum_number_n.value_or(0),
                      quantum_number_nu.value_or(std::numeric_limits<T>::quiet_NaN()), 0,
                      quantum_number_l.value_or(std::numeric_limits<T>::quiet_NaN()), 0,
                      quantum_number_s.value_or(std::numeric_limits<T>::quiet_NaN()), 0,
                      quantum_number_j.value_or(std::numeric_limits<T>::quiet_NaN()), 0);
}

// Explicit instantiations
template class KetAtomCreator<float>;
template class KetAtomCreator<double>;

///////////////////////////////////////////////////////////////////////////////////////
// Test cases
///////////////////////////////////////////////////////////////////////////////////////

#include <doctest/doctest.h>

DOCTEST_TEST_CASE("create a ket for rubidium") {
    auto ket = KetAtomCreator<float>("Rb", 60, 1, 0.5, 0.5).create();
    CHECK(ket.get_species() == "Rb");
    CHECK(ket.get_quantum_number_n() == 60);
    CHECK(ket.get_quantum_number_l() == 1);
    CHECK(ket.get_quantum_number_f() == 0.5);
    CHECK(ket.get_quantum_number_j() == 0.5);
    CHECK(ket.get_quantum_number_m() == 0.5);
    CHECK(ket.get_quantum_number_s() == 0.5);
    CHECK(ket.get_parity() == -1);
}

DOCTEST_TEST_CASE("create a ket for strontium") {
    KetAtomCreator<float> creator("Sr88_mqdt");
    creator.set_quantum_number_nu(60);
    creator.set_quantum_number_l(1);
    creator.set_quantum_number_f(1);
    creator.set_quantum_number_m(0);
    creator.set_quantum_number_s(0);
    auto ket = creator.create();
    CHECK(ket.get_species() == "Sr88_mqdt");
    CHECK(ket.get_quantum_number_nu() == 60);
    CHECK(ket.get_quantum_number_f() == 1);
    CHECK(ket.get_quantum_number_m() == 0);
    CHECK(ket.get_parity() == -1);
}
