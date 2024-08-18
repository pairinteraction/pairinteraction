#include "pintr/ket/KetAtomCreator.hpp"

#include "pintr/database/AtomDescriptionByParameters.hpp"
#include "pintr/database/Database.hpp"
#include "pintr/enums/Parity.hpp"

#include <cmath>

namespace pintr {
template <typename Real>
KetAtomCreator<Real>::KetAtomCreator(std::string species, int n, Real l, Real j, Real m)
    : species(std::move(species)), quantum_number_f(j), quantum_number_m(m), quantum_number_n(n),
      quantum_number_l(l), quantum_number_s(0.5), quantum_number_j(j) {}

template <typename Real>
KetAtomCreator<Real> &KetAtomCreator<Real>::set_species(const std::string &value) {
    species.emplace(value);
    return *this;
}

template <typename Real>
KetAtomCreator<Real> &KetAtomCreator<Real>::set_energy(Real value) {
    energy.emplace(value);
    return *this;
}

template <typename Real>
KetAtomCreator<Real> &KetAtomCreator<Real>::set_quantum_number_f(Real value) {
    if (2 * value != std::rintf(2 * value)) {
        throw std::invalid_argument("Quantum number f must be an integer or half-integer.");
    }
    quantum_number_f.emplace(value);
    return *this;
}

template <typename Real>
KetAtomCreator<Real> &KetAtomCreator<Real>::set_quantum_number_m(Real value) {
    if (2 * value != std::rintf(2 * value)) {
        throw std::invalid_argument("Quantum number m must be an integer or half-integer.");
    }
    quantum_number_m.emplace(value);
    return *this;
}

template <typename Real>
KetAtomCreator<Real> &KetAtomCreator<Real>::set_parity(Parity value) {
    parity = value;
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
std::shared_ptr<const KetAtom<Real>> KetAtomCreator<Real>::create(Database &database) const {

    if (!species.has_value()) {
        throw std::runtime_error("Species not set.");
    }

    AtomDescriptionByParameters<Real> description{
        parity,           energy,           quantum_number_f,
        quantum_number_m, quantum_number_n, quantum_number_nu,
        quantum_number_l, quantum_number_s, quantum_number_j};

    return database.get_ket<Real>(species.value(), description);
}

// Explicit instantiations
template class KetAtomCreator<float>;
template class KetAtomCreator<double>;
} // namespace pintr
