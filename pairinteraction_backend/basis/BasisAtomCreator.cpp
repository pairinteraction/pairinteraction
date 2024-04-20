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
BasisAtomCreator<Scalar> &BasisAtomCreator<Scalar>::add_ket(const ket_t &ket) {
    if (additional_ket_species.has_value() && additional_ket_species.value() != ket.get_species()) {
        throw std::invalid_argument("Species mismatch.");
    }
    additional_ket_species.emplace(ket.get_species());
    additional_ket_ids.emplace_back(ket.get_id());
    return *this;
}

template <typename Scalar>
BasisAtom<Scalar> BasisAtomCreator<Scalar>::create(Database &database) const {

    if (species.has_value() && additional_ket_species.has_value() &&
        species.value() != additional_ket_species.value()) {
        throw std::invalid_argument("Species mismatch.");
    }

    std::string extracted_species;
    if (species.has_value()) {
        extracted_species = species.value();
    } else if (additional_ket_species.has_value()) {
        extracted_species = additional_ket_species.value();
    } else {
        throw std::runtime_error("Species not set.");
    }

    auto kets = database.get_kets<real_t>(
        extracted_species, min_energy, max_energy, min_quantum_number_f, max_quantum_number_f,
        min_quantum_number_m, max_quantum_number_m, parity, min_quantum_number_n,
        max_quantum_number_n, min_quantum_number_nu, max_quantum_number_nu, min_quantum_number_l,
        max_quantum_number_l, min_quantum_number_s, max_quantum_number_s, min_quantum_number_j,
        max_quantum_number_j, additional_ket_ids);

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

#include "ket/KetAtomCreator.hpp"
#include <doctest/doctest.h>
#include <fmt/ostream.h>
#include <spdlog/spdlog.h>

#if FMT_VERSION < 90000
namespace fmt {
template <typename T>
inline auto streamed(T &&v) {
    return std::forward<T>(v);
}
} // namespace fmt
#endif

DOCTEST_TEST_CASE("create a basis for strontium") {
    auto database = Database();
    auto basis = BasisAtomCreator<float>()
                     .set_species("Sr88_singlet")
                     .restrict_quantum_number_n(60, 64)
                     .restrict_quantum_number_l(0, 2)
                     .create(database);
    for (const auto &ket : basis) {
        DOCTEST_CHECK(ket.get_species() == "Sr88_singlet");
        SPDLOG_LOGGER_INFO(spdlog::get("doctest"), "Ket: {}", fmt::streamed(ket));
    }
}

DOCTEST_TEST_CASE("create a basis from kets") {
    auto database = Database();
    auto ket1 = KetAtomCreator<float>("Sr88_singlet", 80, 0, 0, 0).create(database);
    auto ket2 = KetAtomCreator<float>("Sr88_singlet", 81, 0, 0, 0).create(database);
    auto ket3 = KetAtomCreator<float>("Sr88_singlet", 82, 0, 0, 0).create(database);
    auto basis =
        BasisAtomCreator<float>().add_ket(ket1).add_ket(ket2).add_ket(ket3).create(database);
    for (const auto &ket : basis) {
        DOCTEST_CHECK(ket.get_species() == "Sr88_singlet");
        SPDLOG_LOGGER_INFO(spdlog::get("doctest"), "Ket: {}", fmt::streamed(ket));
    }
}
