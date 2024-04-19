#include "basis/BasisAtomCreator.hpp"
#include "basis/BasisAtom.hpp"
#include "database/Database.hpp"

template <typename Scalar>
BasisAtomCreator<Scalar> &BasisAtomCreator<Scalar>::set_species(std::string value) {
    species.emplace(value);
    return *this;
}

template <typename Scalar>
BasisAtomCreator<Scalar> &BasisAtomCreator<Scalar>::restrict_energy(real_t min, real_t max) {
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
BasisAtomCreator<Scalar> &BasisAtomCreator<Scalar>::restrict_quantum_number_nu(real_t min,
                                                                               real_t max) {
    min_quantum_number_nu.emplace(min);
    max_quantum_number_nu.emplace(max);
    return *this;
}

template <typename Scalar>
BasisAtomCreator<Scalar> &BasisAtomCreator<Scalar>::restrict_quantum_number_l(real_t min,
                                                                              real_t max) {
    min_quantum_number_l.emplace(min);
    max_quantum_number_l.emplace(max);
    return *this;
}

template <typename Scalar>
BasisAtomCreator<Scalar> &BasisAtomCreator<Scalar>::restrict_quantum_number_s(real_t min,
                                                                              real_t max) {
    min_quantum_number_s.emplace(min);
    max_quantum_number_s.emplace(max);
    return *this;
}

template <typename Scalar>
BasisAtomCreator<Scalar> &BasisAtomCreator<Scalar>::restrict_quantum_number_j(real_t min,
                                                                              real_t max) {
    min_quantum_number_j.emplace(min);
    max_quantum_number_j.emplace(max);
    return *this;
}

template <typename Scalar>
BasisAtom<Scalar> BasisAtomCreator<Scalar>::create(Database &database) const {

    if (!species.has_value()) {
        throw std::runtime_error("Species not set.");
    }

    auto kets = database.get_kets<real_t>(
        species.value(), min_energy, max_energy, min_quantum_number_f, max_quantum_number_f,
        min_quantum_number_m, max_quantum_number_m, parity, min_quantum_number_n,
        max_quantum_number_n, min_quantum_number_nu, max_quantum_number_nu, min_quantum_number_l,
        max_quantum_number_l, min_quantum_number_s, max_quantum_number_s, min_quantum_number_j,
        max_quantum_number_j);

    return BasisAtom<Scalar>(std::move(kets), database);
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
    auto database = Database();
    auto basis = BasisAtomCreator<float>()
                     .set_species("Rb")
                     .restrict_quantum_number_n(60, 64)
                     .restrict_quantum_number_l(0, 2)
                     .create(database);
    for (const auto &ket : basis) {
        DOCTEST_CHECK(ket.get_species() == "Rb");
    }
}
