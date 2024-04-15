#include "basis/Basis.hpp"

template <typename Derived>
Basis<Derived>::Basis(ketvec_t &&kets) : kets(std::move(kets)), is_standard_basis(true) {
    energies.reserve(kets.size());
    quantum_numbers_f.reserve(kets.size());
    quantum_numbers_m.reserve(kets.size());
    parities.reserve(kets.size());
    labels.reserve(kets.size());
    for (const auto &ket : kets) {
        energies.push_back(ket->get_energy());
        quantum_numbers_f.push_back(ket->get_quantum_number_f());
        quantum_numbers_m.push_back(ket->get_quantum_number_m());
        parities.push_back(ket->get_parity());
        labels.push_back(ket->get_label());
    }
    coefficients = Eigen::SparseMatrix<scalar_t>(kets.size(), kets.size());
    coefficients.setIdentity();
}

template <typename Derived>
const Derived &Basis<Derived>::derived() const {
    return static_cast<const Derived &>(*this);
}

template <typename Derived>
const typename Basis<Derived>::ket_t &Basis<Derived>::get_ket(size_t index) const {
    return *kets[index];
}

template <typename Derived>
typename Basis<Derived>::real_t Basis<Derived>::get_energy(size_t index) const {
    return energies[index];
}

template <typename Derived>
float Basis<Derived>::get_quantum_number_f(size_t index) const {
    return quantum_numbers_f[index];
}

template <typename Derived>
float Basis<Derived>::get_quantum_number_m(size_t index) const {
    return quantum_numbers_m[index];
}

template <typename Derived>
int Basis<Derived>::get_parity(size_t index) const {
    return parities[index];
}

template <typename Derived>
std::string Basis<Derived>::get_label(size_t index) const {
    return labels[index];
}

template <typename Derived>
typename Basis<Derived>::Iterator Basis<Derived>::begin() const {
    return Iterator(*this, 0);
}

template <typename Derived>
typename Basis<Derived>::Iterator Basis<Derived>::end() const {
    return Iterator(*this, kets.size());
}

template <typename Derived>
Basis<Derived>::Iterator::Iterator(const Basis<Derived> &basis, size_t index)
    : basis(basis), index(index) {}

template <typename Derived>
bool Basis<Derived>::Iterator::operator!=(const Iterator &other) const {
    return index != other.index;
}

template <typename Derived>
const typename Basis<Derived>::ket_t &Basis<Derived>::Iterator::operator*() const {
    return basis.get_ket(index);
}

template <typename Derived>
typename Basis<Derived>::Iterator &Basis<Derived>::Iterator::operator++() {
    ++index;
    return *this;
}

template <typename Derived>
std::vector<int> Basis<Derived>::get_sorter_according_to_kets() const {
    throw std::runtime_error("Not implemented");
}

template <typename Derived>
std::vector<int> Basis<Derived>::get_sorter_according_to_energies() const {
    throw std::runtime_error("Not implemented");
}

template <typename Derived>
std::vector<int> Basis<Derived>::get_sorter_according_to_quantum_number_m() const {
    throw std::runtime_error("Not implemented");
}

template <typename Derived>
std::vector<int> Basis<Derived>::get_sorter_according_to_parity() const {
    throw std::runtime_error("Not implemented");
}

template <typename Derived>
void Basis<Derived>::transform(const Eigen::SparseMatrix<scalar_t> &transformator) {
    throw std::runtime_error("Not implemented");
}

template <typename Derived>
void Basis<Derived>::sort(const std::vector<int> &sorter) {
    throw std::runtime_error("Not implemented");
}

template <typename Derived>
void Basis<Derived>::sort_according_to_kets() {
    auto sorter = this->get_sorter_according_to_kets();
    this->sort(sorter);
}

template <typename Derived>
void Basis<Derived>::sort_according_to_energies() {
    auto sorter = this->get_sorter_according_to_energies();
    this->sort(sorter);
}

template <typename Derived>
void Basis<Derived>::sort_according_to_quantum_number_m() {
    auto sorter = this->get_sorter_according_to_quantum_number_m();
    this->sort(sorter);
}

template <typename Derived>
void Basis<Derived>::sort_according_to_parity() {
    auto sorter = this->get_sorter_according_to_parity();
    this->sort(sorter);
}

// Explicit instantiations
#include "basis/BasisAtom.hpp"
#include "basis/BasisClassicalLight.hpp"

template class Basis<BasisAtom<float>>;
template class Basis<BasisAtom<double>>;
template class Basis<BasisAtom<std::complex<float>>>;
template class Basis<BasisAtom<std::complex<double>>>;
template class Basis<BasisClassicalLight<float>>;
template class Basis<BasisClassicalLight<double>>;
template class Basis<BasisClassicalLight<std::complex<float>>>;
template class Basis<BasisClassicalLight<std::complex<double>>>;

///////////////////////////////////////////////////////////////////////////////////////
// Test cases
///////////////////////////////////////////////////////////////////////////////////////

#include <doctest/doctest.h>

#ifndef DOCTEST_CONFIG_DISABLE

// A derived ket class
class KetDerived : public Ket<float> {
public:
    int get_new_property() const { return new_property; }

private:
    friend class KetDerivedCreator;
    KetDerived(float energy, float f, float m, int p, std::string label, int new_property)
        : Ket<float>(energy, f, m, p, label), new_property(new_property) {}
    int new_property;
};

// Classes for creating an instance of the derived ket class
class KetDerivedCreator {
public:
    KetDerivedCreator(float energy, float f, float m, int p, std::string label, int new_property)
        : energy(energy), f(f), m(m), p(p), label(label), new_property(new_property) {}
    KetDerived create() const { return KetDerived(energy, f, m, p, label, new_property); }

private:
    float energy;
    float f;
    float m;
    int p;
    std::string label;
    int new_property;
};

// A derived basis class
class BasisDerived;

template <>
struct internal::BasisTraits<BasisDerived> {
    using scalar_t = float;
    using real_t = float;
    using ket_t = KetDerived;
    using ketvec_t = std::vector<std::shared_ptr<const ket_t>>;
};

class BasisDerived : public Basis<BasisDerived> {
public:
    using Type = BasisDerived;
    using ketvec_t = typename internal::BasisTraits<BasisDerived>::ketvec_t;

private:
    friend class BasisDerivedCreator;
    BasisDerived(ketvec_t &&kets) : Basis<BasisDerived>(std::move(kets)) {}
};

// Classes for creating an instance of the derived basis class
class BasisDerivedCreator {
public:
    BasisDerivedCreator() = default;
    BasisDerived create() const {
        std::vector<std::shared_ptr<const KetDerived>> kets;
        kets.reserve(3);
        kets.push_back(std::make_shared<const KetDerived>(
            KetDerivedCreator(1.0, 0.5, 0.5, 1, "1s", 42).create()));
        kets.push_back(std::make_shared<const KetDerived>(
            KetDerivedCreator(2.0, 0.5, 0.5, 1, "2s", 42).create()));
        kets.push_back(std::make_shared<const KetDerived>(
            KetDerivedCreator(3.0, 0.5, 0.5, 1, "3s", 42).create()));
        return BasisDerived(std::move(kets));
    }
};

#endif // DOCTEST_CONFIG_DISABLE

DOCTEST_TEST_CASE("constructing a class derived from basis") {
    auto basis = BasisDerivedCreator().create();

    // Check that the kets can be iterated over and the new property can be obtained
    for (const auto &ket : basis) {
        DOCTEST_CHECK(ket.get_new_property() == 42);
    }
}
