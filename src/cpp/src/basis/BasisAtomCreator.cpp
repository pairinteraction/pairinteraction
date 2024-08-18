#include "pintr/basis/BasisAtomCreator.hpp"

#include "pintr/database/AtomDescriptionByRanges.hpp"
#include "pintr/database/Database.hpp"
#include "pintr/enums/Parity.hpp"
#include "pintr/ket/KetAtom.hpp"

namespace pintr {
template <typename Scalar>
BasisAtomCreator<Scalar> &BasisAtomCreator<Scalar>::set_species(const std::string &value) {
    species.emplace(value);
    return *this;
}

template <typename Scalar>
BasisAtomCreator<Scalar> &BasisAtomCreator<Scalar>::restrict_energy(real_t min, real_t max) {
    range_energy = {min, max};
    return *this;
}

template <typename Scalar>
BasisAtomCreator<Scalar> &BasisAtomCreator<Scalar>::restrict_quantum_number_f(real_t min,
                                                                              real_t max) {
    if (2 * min != std::rintf(2 * min) || 2 * max != std::rintf(2 * max)) {
        throw std::invalid_argument("Quantum number f must be an integer or half-integer.");
    }
    range_quantum_number_f = {min, max};
    return *this;
}

template <typename Scalar>
BasisAtomCreator<Scalar> &BasisAtomCreator<Scalar>::restrict_quantum_number_m(real_t min,
                                                                              real_t max) {
    if (2 * min != std::rintf(2 * min) || 2 * max != std::rintf(2 * max)) {
        throw std::invalid_argument("Quantum number m must be an integer or half-integer.");
    }
    range_quantum_number_m = {min, max};
    return *this;
}

template <typename Scalar>
BasisAtomCreator<Scalar> &BasisAtomCreator<Scalar>::restrict_parity(Parity value) {
    parity = value;
    return *this;
}

template <typename Scalar>
BasisAtomCreator<Scalar> &BasisAtomCreator<Scalar>::restrict_quantum_number_n(int min, int max) {
    range_quantum_number_n = {min, max};
    return *this;
}

template <typename Scalar>
BasisAtomCreator<Scalar> &BasisAtomCreator<Scalar>::restrict_quantum_number_nu(real_t min,
                                                                               real_t max) {
    range_quantum_number_nu = {min, max};
    return *this;
}

template <typename Scalar>
BasisAtomCreator<Scalar> &BasisAtomCreator<Scalar>::restrict_quantum_number_l(real_t min,
                                                                              real_t max) {
    range_quantum_number_l = {min, max};
    return *this;
}

template <typename Scalar>
BasisAtomCreator<Scalar> &BasisAtomCreator<Scalar>::restrict_quantum_number_s(real_t min,
                                                                              real_t max) {
    range_quantum_number_s = {min, max};
    return *this;
}

template <typename Scalar>
BasisAtomCreator<Scalar> &BasisAtomCreator<Scalar>::restrict_quantum_number_j(real_t min,
                                                                              real_t max) {
    range_quantum_number_j = {min, max};
    return *this;
}

template <typename Scalar>
BasisAtomCreator<Scalar> &BasisAtomCreator<Scalar>::add_ket(std::shared_ptr<const ket_t> ket) {
    if (additional_ket_species.has_value() &&
        additional_ket_species.value() != ket->get_species()) {
        throw std::invalid_argument("Species mismatch.");
    }
    additional_ket_species.emplace(ket->get_species());
    additional_ket_ids.push_back(ket->get_id());
    return *this;
}

template <typename Scalar>
std::shared_ptr<const BasisAtom<Scalar>>
BasisAtomCreator<Scalar>::create(Database &database) const {

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

    AtomDescriptionByRanges<real_t> description{parity,
                                                range_energy,
                                                range_quantum_number_f,
                                                range_quantum_number_m,
                                                range_quantum_number_n,
                                                range_quantum_number_nu,
                                                range_quantum_number_l,
                                                range_quantum_number_s,
                                                range_quantum_number_j};

    return database.get_basis<Scalar>(extracted_species, description, additional_ket_ids);
}

// Explicit instantiations
template class BasisAtomCreator<float>;
template class BasisAtomCreator<double>;
template class BasisAtomCreator<std::complex<float>>;
template class BasisAtomCreator<std::complex<double>>;
} // namespace pintr
