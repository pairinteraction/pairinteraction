#include "ket/KetAtomCreator.hpp"

#include "database/Database.hpp"
#include "ket/KetAtom.hpp"

#include <cmath>
#include <limits>

template <typename Real>
KetAtomCreator<Real>::KetAtomCreator(std::string species, int n, Real l, float j, float m)
    : species(species), quantum_number_f(j), quantum_number_m(m), quantum_number_n(n),
      quantum_number_l(l), quantum_number_s(0.5), quantum_number_j(j) {}

template <typename Real>
KetAtomCreator<Real> &KetAtomCreator<Real>::set_species(std::string value) {
    species.emplace(value);
    return *this;
}

template <typename Real>
KetAtomCreator<Real> &KetAtomCreator<Real>::set_energy(Real value) {
    energy.emplace(value);
    return *this;
}

template <typename Real>
KetAtomCreator<Real> &KetAtomCreator<Real>::set_quantum_number_f(float value) {
    quantum_number_f.emplace(value);
    return *this;
}

template <typename Real>
KetAtomCreator<Real> &KetAtomCreator<Real>::set_quantum_number_m(float value) {
    quantum_number_m.emplace(value);
    return *this;
}

template <typename Real>
KetAtomCreator<Real> &KetAtomCreator<Real>::set_parity(int value) {
    parity.emplace(value);
    return *this;
}

template <typename Real>
KetAtomCreator<Real> &KetAtomCreator<Real>::set_quantum_number_n(int value) {
    quantum_number_n.emplace(value);
    return *this;
}

template <typename Real>
KetAtomCreator<Real> &KetAtomCreator<Real>::set_quantum_number_nu(Real value) {
    quantum_number_nu.emplace(value);
    return *this;
}

template <typename Real>
KetAtomCreator<Real> &KetAtomCreator<Real>::set_quantum_number_l(Real value) {
    quantum_number_l.emplace(value);
    return *this;
}

template <typename Real>
KetAtomCreator<Real> &KetAtomCreator<Real>::set_quantum_number_s(Real value) {
    quantum_number_s.emplace(value);
    return *this;
}

template <typename Real>
KetAtomCreator<Real> &KetAtomCreator<Real>::set_quantum_number_j(Real value) {
    quantum_number_j.emplace(value);
    return *this;
}

template <typename Real>
KetAtom<Real> KetAtomCreator<Real>::create(Database &database) const {

    if (!species.has_value()) {
        throw std::runtime_error("Species not set.");
    }

    auto ket = database.get_ket<Real>(species.value(), energy, quantum_number_f, quantum_number_m,
                                      parity, quantum_number_n, quantum_number_nu, quantum_number_l,
                                      quantum_number_s, quantum_number_j);

    return ket;
}

// Explicit instantiations
template class KetAtomCreator<float>;
template class KetAtomCreator<double>;

///////////////////////////////////////////////////////////////////////////////////////
// Test cases
///////////////////////////////////////////////////////////////////////////////////////

#include "utils/streamed.hpp"
#include <doctest/doctest.h>
#include <spdlog/spdlog.h>

// DOCTEST_TEST_CASE("create a ket for rubidium") {
//     auto database = Database();
//     auto ket = KetAtomCreator<float>("Rb", 60, 1, 0.5, 0.5).create(database);
//     DOCTEST_CHECK(ket.get_species() == "Rb");
//     DOCTEST_CHECK(ket.get_quantum_number_n() == 60);
//     DOCTEST_CHECK(ket.get_quantum_number_l() == 1);
//     DOCTEST_CHECK(ket.get_quantum_number_f() == 0.5);
//     DOCTEST_CHECK(ket.get_quantum_number_j() == 0.5);
//     DOCTEST_CHECK(ket.get_quantum_number_m() == 0.5);
//     DOCTEST_CHECK(ket.get_quantum_number_s() == 0.5);
//     DOCTEST_CHECK(ket.get_parity() == -1);
//     SPDLOG_LOGGER_INFO(spdlog::get("doctest"), "Ket: {}", fmt::streamed(ket));
// }

DOCTEST_TEST_CASE("create a ket for strontium") {
    auto database = Database();
    auto ket = KetAtomCreator<float>()
                   .set_species("Sr88_singlet")
                   .set_quantum_number_n(60)
                   .set_quantum_number_l(1)
                   .set_quantum_number_f(1)
                   .set_quantum_number_m(0)
                   .set_quantum_number_s(0)
                   .create(database);
    DOCTEST_CHECK(ket.get_species() == "Sr88_singlet");
    DOCTEST_CHECK(ket.get_quantum_number_n() == 60);
    DOCTEST_CHECK(ket.get_quantum_number_f() == 1);
    DOCTEST_CHECK(ket.get_quantum_number_m() == 0);
    DOCTEST_CHECK(ket.get_parity() == -1);
    SPDLOG_LOGGER_INFO(spdlog::get("doctest"), "Ket: {}", fmt::streamed(ket));
}
