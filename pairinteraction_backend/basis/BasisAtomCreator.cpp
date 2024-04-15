#include "basis/BasisAtomCreator.hpp"

template <typename Scalar>
BasisAtomCreator<Scalar>::BasisAtomCreator(std::string species) : species(species) {}

template <typename Scalar>
BasisAtomCreator<Scalar> &BasisAtomCreator<Scalar>::restrict_energy(Scalar min, Scalar max) {
    min_energy.emplace(min);
    max_energy.emplace(max);
    return *this;
}

template <typename Scalar>
BasisAtomCreator<Scalar> &BasisAtomCreator<Scalar>::restrict_quantum_number_f(float min,
                                                                              float max) {
    min_quantum_number_f.emplace(min);
    max_quantum_number_f.emplace(max);
    return *this;
}

template <typename Scalar>
BasisAtomCreator<Scalar> &BasisAtomCreator<Scalar>::restrict_quantum_number_m(float min,
                                                                              float max) {
    min_quantum_number_m.emplace(min);
    max_quantum_number_m.emplace(max);
    return *this;
}

template <typename Scalar>
BasisAtomCreator<Scalar> &BasisAtomCreator<Scalar>::restrict_parity(int parity) {
    this->parity.emplace(parity);
    return *this;
}

template <typename Scalar>
BasisAtomCreator<Scalar> &BasisAtomCreator<Scalar>::restrict_quantum_number_n(int min, int max) {
    min_quantum_number_n.emplace(min);
    max_quantum_number_n.emplace(max);
    return *this;
}

template <typename Scalar>
BasisAtomCreator<Scalar> &BasisAtomCreator<Scalar>::restrict_quantum_number_nu(Scalar min,
                                                                               Scalar max) {
    min_quantum_number_nu.emplace(min);
    max_quantum_number_nu.emplace(max);
    return *this;
}

template <typename Scalar>
BasisAtomCreator<Scalar> &BasisAtomCreator<Scalar>::restrict_quantum_number_l(Scalar min,
                                                                              Scalar max) {
    min_quantum_number_l.emplace(min);
    max_quantum_number_l.emplace(max);
    return *this;
}

template <typename Scalar>
BasisAtomCreator<Scalar> &BasisAtomCreator<Scalar>::restrict_quantum_number_s(Scalar min,
                                                                              Scalar max) {
    min_quantum_number_s.emplace(min);
    max_quantum_number_s.emplace(max);
    return *this;
}

template <typename Scalar>
BasisAtomCreator<Scalar> &BasisAtomCreator<Scalar>::restrict_quantum_number_j(Scalar min,
                                                                              Scalar max) {
    min_quantum_number_j.emplace(min);
    max_quantum_number_j.emplace(max);
    return *this;
}

template <typename Scalar>
BasisAtom<Scalar> BasisAtomCreator<Scalar>::create() const {

    // TODO perform database request

    std::vector<std::shared_ptr<const ket_t>> kets;
    kets.reserve(2);
    kets.push_back(
        std::make_shared<const ket_t>(KetAtomCreator<real_t>(species, 60, 1, 0.5, -0.5).create()));
    kets.push_back(
        std::make_shared<const ket_t>(KetAtomCreator<real_t>(species, 60, 1, 0.5, 0.5).create()));

    return BasisAtom<Scalar>(std::move(kets));
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
    auto basis = BasisAtomCreator<float>("Rb")
                     .restrict_quantum_number_n(60, 64)
                     .restrict_quantum_number_l(0, 2)
                     .create();
    for (const auto &ket : basis) {
        DOCTEST_CHECK(ket.get_species() == "Rb");
    }
}
