#include "Basis.hpp"

template <typename T, bool is_complex>
Basis<T, is_complex>::Basis() : is_assembled(false), is_standard_basis(true) {}

template <typename T, bool is_complex>
const Ket<T> &Basis<T, is_complex>::get_ket(size_t index) {
    this->ensure_assembled();
    return *kets[index];
}

template <typename T, bool is_complex>
T Basis<T, is_complex>::get_energy(size_t index) {
    this->ensure_assembled();
    return energies[index];
}

template <typename T, bool is_complex>
float Basis<T, is_complex>::get_quantum_number_f(size_t index) {
    this->ensure_assembled();
    return quantum_numbers_f[index];
}

template <typename T, bool is_complex>
float Basis<T, is_complex>::get_quantum_number_m(size_t index) {
    this->ensure_assembled();
    return quantum_numbers_m[index];
}

template <typename T, bool is_complex>
int Basis<T, is_complex>::get_parity(size_t index) {
    this->ensure_assembled();
    return parities[index];
}

template <typename T, bool is_complex>
std::string Basis<T, is_complex>::get_label(size_t index) {
    this->ensure_assembled();
    return labels[index];
}

template <typename T, bool is_complex>
void Basis<T, is_complex>::ensure_assembled() {
    if (!is_assembled) {
        this->ensure_assembled_kets();

        energies.reserve(kets.size());
        quantum_numbers_f.reserve(kets.size());
        quantum_numbers_m.reserve(kets.size());
        parities.reserve(kets.size());
        labels.reserve(kets.size());
        for (auto &ket : kets) {
            energies.push_back(ket->get_energy());
            quantum_numbers_f.push_back(ket->get_quantum_number_f());
            quantum_numbers_m.push_back(ket->get_quantum_number_m());
            parities.push_back(ket->get_parity());
            labels.push_back(ket->get_label());
        }
        coefficients = Eigen::SparseMatrix<scalar_t>(kets.size(), kets.size());
        coefficients.setIdentity();
    }

    is_assembled = true;
}

template <typename T, bool is_complex>
void Basis<T, is_complex>::ensure_not_assembled() const {
    if (is_assembled) {
        throw std::runtime_error("Basis is already assembled");
    }
}

template <typename T, bool is_complex>
void Basis<T, is_complex>::ensure_standard_basis() const {
    if (!is_standard_basis) {
        throw std::runtime_error("Basis is not a standard basis");
    }
}

template <typename T, bool is_complex>
typename Basis<T, is_complex>::Iterator Basis<T, is_complex>::begin() {
    this->ensure_assembled();
    return Iterator(*this, 0);
}

template <typename T, bool is_complex>
typename Basis<T, is_complex>::Iterator Basis<T, is_complex>::end() {
    this->ensure_assembled();
    return Iterator(*this, kets.size());
}

template <typename T, bool is_complex>
Basis<T, is_complex>::Iterator::Iterator(const Basis<T, is_complex> &basis, size_t index)
    : basis(basis), index(index) {}

template <typename T, bool is_complex>
bool Basis<T, is_complex>::Iterator::operator!=(const Iterator &other) const {
    return index != other.index;
}

template <typename T, bool is_complex>
const Ket<T> &Basis<T, is_complex>::Iterator::operator*() const {
    return *basis.kets[index];
}

template <typename T, bool is_complex>
typename Basis<T, is_complex>::Iterator &Basis<T, is_complex>::Iterator::operator++() {
    ++index;
    return *this;
}

// Explicit instantiations
template class Basis<float, false>;
template class Basis<double, false>;
template class Basis<float, true>;
template class Basis<double, true>;

///////////////////////////////////////////////////////////////////////////////////////
// Test cases
///////////////////////////////////////////////////////////////////////////////////////

#include <doctest/doctest.h>

DOCTEST_TEST_CASE("constructing a class derived from basis") {

    class KetDerived : public Ket<float> {
    public:
        KetDerived(float energy, float f, float m, int p, std::string label, int new_property)
            : Ket<float>(energy, f, m, p, label), new_property(new_property) {}
        int get_new_property() const { return new_property; }

    private:
        void ensure_assembled_ket() override {}
        int new_property;
    };

    class BasisDerived : public Basis<float, false> {
    public:
        void ensure_assembled_kets() override {
            this->kets.reserve(3);
            this->kets.push_back(std::make_shared<KetDerived>(1.0, 0.5, 0.5, 1, "1s", 42));
            this->kets.push_back(std::make_shared<KetDerived>(2.0, 0.5, 0.5, 1, "2s", 42));
            this->kets.push_back(std::make_shared<KetDerived>(3.0, 0.5, 0.5, 1, "3s", 42));
        }
    };

    BasisDerived basis;

    // Check that the kets can be iterated over and the new property can be obtained
    for (const auto &ket : basis) {
        DOCTEST_CHECK(dynamic_cast<const KetDerived &>(ket).get_new_property() == 42);
    }
}
