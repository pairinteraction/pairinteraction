#include "Basis.hpp"

template <typename T, bool is_complex>
Basis<T, is_complex>::Basis() : is_assembled(false), is_standard_basis(true), num_coefficients(0) {}

template <typename T, bool is_complex>
T Basis<T, is_complex>::get_energy(size_t index) const {
    if (!is_assembled) {
        throw std::runtime_error("The basis has not been assembled");
    }
    return energies[index];
}

template <typename T, bool is_complex>
float Basis<T, is_complex>::get_quantum_number_f(size_t index) const {
    if (!is_assembled) {
        throw std::runtime_error("The basis has not been assembled");
    }
    return quantum_numbers_f[index];
}

template <typename T, bool is_complex>
float Basis<T, is_complex>::get_quantum_number_m(size_t index) const {
    if (!is_assembled) {
        throw std::runtime_error("The basis has not been assembled");
    }
    return quantum_numbers_m[index];
}

template <typename T, bool is_complex>
int Basis<T, is_complex>::get_parity(size_t index) const {
    if (!is_assembled) {
        throw std::runtime_error("The basis has not been assembled");
    }
    return parities[index];
}

template <typename T, bool is_complex>
std::string Basis<T, is_complex>::get_label(size_t index) const {
    if (!is_assembled) {
        throw std::runtime_error("The basis has not been assembled");
    }
    return labels[index];
}

template <typename T, bool is_complex>
void Basis<T, is_complex>::reserve_coefficients(int n) {
    if (is_assembled) {
        throw std::runtime_error("Cannot reserve coefficients after the basis has been assembled");
    }

    num_coefficients = n;
    energies.reserve(n);
    quantum_numbers_f.reserve(n);
    quantum_numbers_m.reserve(n);
    parities.reserve(n);
    labels.reserve(n);
}

template <typename T, bool is_complex>
void Basis<T, is_complex>::add_coefficients(T energy, float f, float m, int p, std::string label) {
    if (is_assembled) {
        throw std::runtime_error("Cannot add coefficients after the basis has assembled");
    }
    if (energies.size() >= num_coefficients) {
        throw std::runtime_error("Cannot add more coefficients than reserved");
    }

    energies.push_back(energy);
    quantum_numbers_f.push_back(f);
    quantum_numbers_m.push_back(m);
    parities.push_back(p);
    labels.push_back(label);
}

template <typename T, bool is_complex>
void Basis<T, is_complex>::assemble_coefficients() {
    if (is_assembled) {
        throw std::runtime_error("The basis has already been assembled");
    }
    if (energies.size() != num_coefficients) {
        throw std::runtime_error("Not all coefficients have been added");
    }

    coefficients = Eigen::SparseMatrix<scalar_t>(num_coefficients, num_coefficients);
    coefficients.setIdentity();
    is_assembled = true;
}

template <typename T, bool is_complex>
typename Basis<T, is_complex>::Iterator Basis<T, is_complex>::begin() const {
    return Iterator(*this, 0);
}

template <typename T, bool is_complex>
typename Basis<T, is_complex>::Iterator Basis<T, is_complex>::end() const {
    return Iterator(*this, num_coefficients);
}

template <typename T, bool is_complex>
Basis<T, is_complex>::Iterator::Iterator(const Basis<T, is_complex> &basis, size_t index)
    : basis(basis), index(index) {}

template <typename T, bool is_complex>
bool Basis<T, is_complex>::Iterator::operator!=(const Iterator &other) const {
    return index != other.index;
}

template <typename T, bool is_complex>
State<T> Basis<T, is_complex>::Iterator::operator*() const {
    return basis.get_state(index);
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
    class StateDerived : public State<float> {
    public:
        StateDerived(float energy, float f, float m, int p, std::string label)
            : State<float>(energy, f, m, p, label) {}
    };

    class BasisDerived : public Basis<float, false> {
    public:
        BasisDerived() : Basis() {
            DOCTEST_CHECK_NOTHROW(this->reserve_coefficients(3));
            DOCTEST_CHECK_NOTHROW(this->add_coefficients(1.0, 0.5, 0.5, 1, "1s"));
            DOCTEST_CHECK_NOTHROW(this->add_coefficients(2.0, 0.5, 0.5, 1, "2s"));
            DOCTEST_CHECK_NOTHROW(this->add_coefficients(3.0, 0.5, 0.5, 1, "3s"));
            DOCTEST_CHECK_NOTHROW(this->assemble_coefficients());

            DOCTEST_CHECK_THROWS_AS(this->reserve_coefficients(3), std::runtime_error);
            DOCTEST_CHECK_THROWS_AS(this->add_coefficients(4.0, 0.5, 0.5, 1, "4s"),
                                    std::runtime_error);
            DOCTEST_CHECK_THROWS_AS(this->assemble_coefficients(), std::runtime_error);
        }
        State<float> get_state(size_t index) const {
            // TODO how to use StateDerived as return type? !!!!!!!!!!!!!!!!!
            return State<float>(this->get_energy(index), this->get_quantum_number_f(index),
                                this->get_quantum_number_m(index), this->get_parity(index),
                                this->get_label(index));
        }
    };

    BasisDerived basis;

    for (const auto &state : basis) {
        DOCTEST_CHECK(state.get_parity() == 1);
    }
}
