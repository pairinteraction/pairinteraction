#include "basis/Basis.hpp"

template <typename T>
Basis<T>::Basis(std::vector<std::shared_ptr<const Ket<real_t<T>>>> &&kets)
    : kets(std::move(kets)), is_standard_basis(true) {
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
    coefficients = Eigen::SparseMatrix<T>(kets.size(), kets.size());
    coefficients.setIdentity();
}

template <typename T>
const Ket<real_t<T>> &Basis<T>::get_ket(size_t index) const {
    return *kets[index];
}

template <typename T>
real_t<T> Basis<T>::get_energy(size_t index) const {
    return energies[index];
}

template <typename T>
float Basis<T>::get_quantum_number_f(size_t index) const {
    return quantum_numbers_f[index];
}

template <typename T>
float Basis<T>::get_quantum_number_m(size_t index) const {
    return quantum_numbers_m[index];
}

template <typename T>
int Basis<T>::get_parity(size_t index) const {
    return parities[index];
}

template <typename T>
std::string Basis<T>::get_label(size_t index) const {
    return labels[index];
}

template <typename T>
typename Basis<T>::Iterator Basis<T>::begin() const {
    return Iterator(*this, 0);
}

template <typename T>
typename Basis<T>::Iterator Basis<T>::end() const {
    return Iterator(*this, kets.size());
}

template <typename T>
Basis<T>::Iterator::Iterator(const Basis<T> &basis, size_t index) : basis(basis), index(index) {}

template <typename T>
bool Basis<T>::Iterator::operator!=(const Iterator &other) const {
    return index != other.index;
}

template <typename T>
const Ket<real_t<T>> &Basis<T>::Iterator::operator*() const {
    return basis.get_ket(index);
}

template <typename T>
typename Basis<T>::Iterator &Basis<T>::Iterator::operator++() {
    ++index;
    return *this;
}

// Explicit instantiations
template class Basis<float>;
template class Basis<double>;
template class Basis<std::complex<float>>;
template class Basis<std::complex<double>>;

///////////////////////////////////////////////////////////////////////////////////////
// Test cases
///////////////////////////////////////////////////////////////////////////////////////

#include <doctest/doctest.h>

DOCTEST_TEST_CASE("constructing a class derived from basis") {

    class KetDerived : public Ket<float> {
    public:
        int get_new_property() const { return new_property; }

    private:
        friend class KetDerivedCreator;
        KetDerived(float energy, float f, float m, int p, std::string label, int new_property)
            : Ket<float>(energy, f, m, p, label), new_property(new_property) {}
        int new_property;
    };

    class KetDerivedCreator {
    public:
        KetDerivedCreator(float energy, float f, float m, int p, std::string label,
                          int new_property)
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

    class BasisDerived : public Basis<float> {
    private:
        friend class BasisDerivedCreator;
        BasisDerived(std::vector<std::shared_ptr<const Ket<float>>> &&kets)
            : Basis<float>(std::move(kets)) {}
    };

    class BasisDerivedCreator {
    public:
        BasisDerivedCreator() = default;
        BasisDerived create() const {
            std::vector<std::shared_ptr<const Ket<float>>> kets;
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

    auto basis = BasisDerivedCreator().create();

    // Check that the kets can be iterated over and the new property can be obtained
    for (const auto &ket : basis) {
        DOCTEST_CHECK(dynamic_cast<const KetDerived &>(ket).get_new_property() == 42);
    }
}
