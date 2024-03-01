#include "basis/Basis.hpp"
#include "basis/BasisAtom.hpp"

template <typename Derived>
Basis<Derived>::Basis() : is_standard_basis(true) {
    energies.reserve(derived().get_kets().size());
    quantum_numbers_f.reserve(derived().get_kets().size());
    quantum_numbers_m.reserve(derived().get_kets().size());
    parities.reserve(derived().get_kets().size());
    labels.reserve(derived().get_kets().size());
    for (const auto &ket : derived().get_kets()) {
        energies.push_back(ket->get_energy());
        quantum_numbers_f.push_back(ket->get_quantum_number_f());
        quantum_numbers_m.push_back(ket->get_quantum_number_m());
        parities.push_back(ket->get_parity());
        labels.push_back(ket->get_label());
    }
    coefficients =
        Eigen::SparseMatrix<Scalar>(derived().get_kets().size(), derived().get_kets().size());
    coefficients.setIdentity();
}

template <typename Derived>
const Derived &Basis<Derived>::derived() const {
    return static_cast<const Derived &>(*this);
}

template <typename Derived>
const auto &Basis<Derived>::get_ket(size_t index) const {
    return derived().get_ket(index); // TODO
}

template <typename Derived>
typename Basis<Derived>::Real Basis<Derived>::get_energy(size_t index) const {
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
    return Iterator(*this, derived().get_kets().size());
}

template <typename Derived>
Basis<Derived>::Iterator::Iterator(const Basis<Derived> &basis, size_t index)
    : basis(basis), index(index) {}

template <typename Derived>
bool Basis<Derived>::Iterator::operator!=(const Iterator &other) const {
    return index != other.index;
}

template <typename Derived>
const Basis<Derived>::KetType &Basis<Derived>::Iterator::operator*() const {
    return basis.get_ket(index);
}

template <typename Derived>
typename Basis<Derived>::Iterator &Basis<Derived>::Iterator::operator++() {
    ++index;
    return *this;
}

// Explicit instantiations
template class Basis<BasisAtom<float>>;
template class Basis<BasisAtom<double>>;
template class Basis<BasisAtom<std::complex<float>>>;
template class Basis<BasisAtom<std::complex<double>>>;

///////////////////////////////////////////////////////////////////////////////////////
// Test cases
///////////////////////////////////////////////////////////////////////////////////////

// #include <doctest/doctest.h>

// DOCTEST_TEST_CASE("constructing a class derived from basis") {

//     class KetDerived : public Ket<float> {
//     public:
//         int get_new_property() const { return new_property; }

//     private:
//         friend class KetDerivedCreator;
//         KetDerived(float energy, float f, float m, int p, std::string label, int new_property)
//             : Ket<float>(energy, f, m, p, label), new_property(new_property) {}
//         int new_property;
//     };

//     class KetDerivedCreator {
//     public:
//         KetDerivedCreator(float energy, float f, float m, int p, std::string label,
//                           int new_property)
//             : energy(energy), f(f), m(m), p(p), label(label), new_property(new_property) {}
//         KetDerived create() const { return KetDerived(energy, f, m, p, label, new_property); }

//     private:
//         float energy;
//         float f;
//         float m;
//         int p;
//         std::string label;
//         int new_property;
//     };

//     class BasisDerived : public Basis<float> {
//     private:
//         friend class BasisDerivedCreator;
//         BasisDerived(std::vector<std::shared_ptr<const Ket<float>>> &&kets)
//             : Basis<float>(std::move(kets)) {}
//     };

//     class BasisDerivedCreator {
//     public:
//         BasisDerivedCreator() = default;
//         BasisDerived create() const {
//             std::vector<std::shared_ptr<const Ket<float>>> kets;
//             kets.reserve(3);
//             kets.push_back(std::make_shared<const KetDerived>(
//                 KetDerivedCreator(1.0, 0.5, 0.5, 1, "1s", 42).create()));
//             kets.push_back(std::make_shared<const KetDerived>(
//                 KetDerivedCreator(2.0, 0.5, 0.5, 1, "2s", 42).create()));
//             kets.push_back(std::make_shared<const KetDerived>(
//                 KetDerivedCreator(3.0, 0.5, 0.5, 1, "3s", 42).create()));
//             return BasisDerived(std::move(kets));
//         }
//     };

//     auto basis = BasisDerivedCreator().create();

//     // Check that the kets can be iterated over and the new property can be obtained
//     for (const auto &ket : basis) {
//         DOCTEST_CHECK(dynamic_cast<const KetDerived &>(ket).get_new_property() == 42);
//     }
// }
