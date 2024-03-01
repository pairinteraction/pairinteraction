#include "basis/BasisAtomCreator.hpp" // TODO all includes needed?
#include "basis/Basis.hpp"
#include "basis/BasisAtom.hpp"
#include "ket/Ket.hpp"
#include "ket/KetAtom.hpp"
#include "ket/KetAtomCreator.hpp"

template <typename T>
BasisAtomCreator<T>::BasisAtomCreator(std::string species) : species(species) {}

template <typename T>
void BasisAtomCreator<T>::restrict_energy(T min, T max) {
    min_energy.emplace(min);
    max_energy.emplace(max);
}

template <typename T>
void BasisAtomCreator<T>::restrict_quantum_number_f(float min, float max) {
    min_quantum_number_f.emplace(min);
    max_quantum_number_f.emplace(max);
}

template <typename T>
void BasisAtomCreator<T>::restrict_quantum_number_m(float min, float max) {
    min_quantum_number_m.emplace(min);
    max_quantum_number_m.emplace(max);
}

template <typename T>
void BasisAtomCreator<T>::restrict_parity(int parity) {
    this->parity.emplace(parity);
}

template <typename T>
void BasisAtomCreator<T>::restrict_quantum_number_n(int min, int max) {
    min_quantum_number_n.emplace(min);
    max_quantum_number_n.emplace(max);
}

template <typename T>
void BasisAtomCreator<T>::restrict_quantum_number_nu(T min, T max) {
    min_quantum_number_nu.emplace(min);
    max_quantum_number_nu.emplace(max);
}

template <typename T>
void BasisAtomCreator<T>::restrict_quantum_number_l(T min, T max) {
    min_quantum_number_l.emplace(min);
    max_quantum_number_l.emplace(max);
}

template <typename T>
void BasisAtomCreator<T>::restrict_quantum_number_s(T min, T max) {
    min_quantum_number_s.emplace(min);
    max_quantum_number_s.emplace(max);
}

template <typename T>
void BasisAtomCreator<T>::restrict_quantum_number_j(T min, T max) {
    min_quantum_number_j.emplace(min);
    max_quantum_number_j.emplace(max);
}

template <typename T>
BasisAtom<T> BasisAtomCreator<T>::create() const {

    // TODO perform database request

    std::vector<std::shared_ptr<const KetAtom<real_t<T>>>> kets;
    kets.reserve(2);
    kets.push_back(std::make_shared<const KetAtom<real_t<T>>>(
        KetAtomCreator<real_t<T>>(species, 60, 1, 0.5, -0.5).create()));
    kets.push_back(std::make_shared<const KetAtom<real_t<T>>>(
        KetAtomCreator<real_t<T>>(species, 60, 1, 0.5, 0.5).create()));

    return BasisAtom<T>(std::move(kets));
}

// Explicit instantiations
template class BasisAtomCreator<float>;
template class BasisAtomCreator<double>;
template class BasisAtomCreator<std::complex<float>>;
template class BasisAtomCreator<std::complex<double>>;

///////////////////////////////////////////////////////////////////////////////////////
// Test cases
///////////////////////////////////////////////////////////////////////////////////////

#include <doctest/doctest.h>

DOCTEST_TEST_CASE("create a basis for rubidium") {
    auto basis = BasisAtomCreator<float>("Rb").create();
    for (const auto &ket : basis) {
        CHECK(ket.get_species() == "Rb");
    }
}
