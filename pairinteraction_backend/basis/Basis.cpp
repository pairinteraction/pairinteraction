#include "basis/Basis.hpp"
#include "utils/Euler.hpp"

#include <numeric>

template <typename Derived>
Basis<Derived>::Basis(ketvec_t &&kets) : kets(std::move(kets)), is_standard_basis(true) {
    energy_of_states.reserve(kets.size());
    for (const auto &ket : kets) {
        energy_of_states.push_back(ket->get_energy());
    }
    coefficients = Eigen::SparseMatrix<scalar_t>(kets.size(), kets.size());
    coefficients.setIdentity();
}

template <typename Derived>
const Derived &Basis<Derived>::derived() const {
    return static_cast<const Derived &>(*this);
}

template <typename Derived>
const typename Basis<Derived>::ket_t &Basis<Derived>::get_ket(size_t index_ket) const {
    return *kets[index_ket];
}

template <typename Derived>
typename Basis<Derived>::real_t Basis<Derived>::get_energy(size_t index_state) const {
    return energy_of_states[index_state];
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
Eigen::SparseMatrix<typename Basis<Derived>::scalar_t>
Basis<Derived>::get_rotator(real_t alpha, real_t beta, real_t gamma) const {
    throw std::runtime_error("Not implemented");
}

template <typename Derived>
Eigen::SparseMatrix<typename Basis<Derived>::scalar_t>
Basis<Derived>::get_rotator(std::array<real_t, 3> to_z_axis,
                            std::array<real_t, 3> to_y_axis) const {
    auto euler_angles = Euler::get_euler_angles(to_z_axis, to_y_axis);
    return this->get_rotator(euler_angles[0], euler_angles[1], euler_angles[2]);
}

template <typename Derived>
std::vector<int> Basis<Derived>::get_sorter_according_to_kets() const {
    std::vector<float> index_ket_of_states;
    throw std::runtime_error("Not implemented");

    std::vector<int> sorter(index_ket_of_states.size());
    std::iota(sorter.begin(), sorter.end(), 0);
    std::stable_sort(sorter.begin(), sorter.end(),
                     [&](int i, int j) { return index_ket_of_states[i] < index_ket_of_states[j]; });
    return sorter;
}

template <typename Derived>
std::vector<int> Basis<Derived>::get_sorter_according_to_energies() const {
    std::vector<int> sorter(energy_of_states.size());
    std::iota(sorter.begin(), sorter.end(), 0);
    std::stable_sort(sorter.begin(), sorter.end(),
                     [&](int i, int j) { return energy_of_states[i] < energy_of_states[j]; });
    return sorter;
}

template <typename Derived>
std::vector<int> Basis<Derived>::get_sorter_according_to_quantum_number_m() const {
    std::vector<float> quantum_number_m_of_states;
    throw std::runtime_error("Not implemented");

    std::vector<int> sorter(quantum_number_m_of_states.size());
    std::iota(sorter.begin(), sorter.end(), 0);
    std::stable_sort(sorter.begin(), sorter.end(), [&](int i, int j) {
        return quantum_number_m_of_states[i] < quantum_number_m_of_states[j];
    });
    return sorter;
}

template <typename Derived>
std::vector<int> Basis<Derived>::get_sorter_according_to_parity() const {
    std::vector<int> parity_of_states;
    throw std::runtime_error("Not implemented");

    std::vector<int> sorter(parity_of_states.size());
    std::iota(sorter.begin(), sorter.end(), 0);
    std::stable_sort(sorter.begin(), sorter.end(),
                     [&](int i, int j) { return parity_of_states[i] < parity_of_states[j]; });
    return sorter;
}

template <typename Derived>
void Basis<Derived>::transform(const Eigen::SparseMatrix<scalar_t> &transformator) {
    coefficients = coefficients * transformator;
}

template <typename Derived>
void Basis<Derived>::rotate(real_t alpha, real_t beta, real_t gamma) {
    auto rotator = this->get_rotator(alpha, beta, gamma);
    this->transform(rotator);
}

template <typename Derived>
void Basis<Derived>::rotate(std::array<real_t, 3> to_z_axis, std::array<real_t, 3> to_y_axis) {
    auto rotator = this->get_rotator(to_z_axis, to_y_axis);
    this->transform(rotator);
}

template <typename Derived>
void Basis<Derived>::sort(const std::vector<int> &sorter) {
    {
        auto tmp(energy_of_states);
        for (int i = 0; i < sorter.size(); ++i) {
            energy_of_states[i] = tmp[sorter[i]];
        }
    }
    Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> perm;
    perm.indices() = Eigen::Map<const Eigen::VectorXi>(sorter.data(), sorter.size());
    coefficients = coefficients * perm;
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

template <typename Derived>
void Basis<Derived>::set_eigen_basis(Eigen::SparseMatrix<scalar_t> evecs,
                                     std::vector<real_t> evals) {
    coefficients = evecs;
    energy_of_states = evals;
    is_standard_basis = false;
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
struct Traits::BasisTraits<BasisDerived> {
    using scalar_t = float;
    using real_t = float;
    using ket_t = KetDerived;
    using ketvec_t = std::vector<std::shared_ptr<const ket_t>>;
};

class BasisDerived : public Basis<BasisDerived> {
public:
    using Type = BasisDerived;
    using ketvec_t = typename Traits::BasisTraits<BasisDerived>::ketvec_t;

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

DOCTEST_TEST_CASE("constructing a class derived from basis") {
    auto basis = BasisDerivedCreator().create();

    // Check that the kets can be iterated over and the new property can be obtained
    for (const auto &ket : basis) {
        DOCTEST_CHECK(ket.get_new_property() == 42);
    }
}

#endif // DOCTEST_CONFIG_DISABLE
