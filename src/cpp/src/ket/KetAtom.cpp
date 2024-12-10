#include "pairinteraction/ket/KetAtom.hpp"

#include "pairinteraction/enums/Parity.hpp"
#include "pairinteraction/utils/hash.hpp"

#include <array>
#include <cctype>
#include <cmath>
#include <fmt/core.h>
#include <string>
#include <string_view>
#include <vector>

namespace pairinteraction {

constexpr std::array<std::string_view, 6> quantum_number_l_labels = {"S", "P", "D", "F", "G", "H"};

template <typename Real>
KetAtom<Real>::KetAtom(Private /*unused*/, Real energy, Real f, Real m, Parity p,
                       std::string species, int n, Real nu_exp, Real nu_std, Real l_exp, Real l_std,
                       Real s_exp, Real s_std, Real j_exp, Real j_std, Database &database,
                       size_t id_in_database)
    : Ket<Real>(energy, f, m, p), species(std::move(species)), quantum_number_n(n),
      quantum_number_nu_exp(nu_exp), quantum_number_nu_std(nu_std), quantum_number_l_exp(l_exp),
      quantum_number_l_std(l_std), quantum_number_s_exp(s_exp), quantum_number_s_std(s_std),
      quantum_number_j_exp(j_exp), quantum_number_j_std(j_std), database(database),
      id_in_database(id_in_database) {}

template <typename Real>
Database &KetAtom<Real>::get_database() const {
    return database;
}

template <typename Real>
size_t KetAtom<Real>::get_id_in_database() const {
    return id_in_database;
}

template <typename Real>
std::string KetAtom<Real>::get_label() const {
    // For the moment, the use of multiple channels is indicated by a negative value of n. When we
    // have updated the database, we must adapt this.
    bool is_single_channel = quantum_number_n > 0;

    size_t pos = species.find('_');
    std::string label = (pos != std::string::npos) ? species.substr(0, pos) : species;
    label[0] = static_cast<char>(std::toupper(label[0]));

    if (is_single_channel) {
        if (quantum_number_s_exp == 0) {
            label += "_singlet";
        } else if (quantum_number_s_exp == 1) {
            label += "_triplet";
        } else if (quantum_number_s_exp != 0.5) {
            throw std::runtime_error(
                "Invalid value for quantum number s in the single-channel description.");
        }
    }

    label += ":";

    if (is_single_channel) {
        label += fmt::format("{:d},", quantum_number_n);
        if (quantum_number_l_exp == std::rintf(quantum_number_l_exp) &&
            quantum_number_l_exp < quantum_number_l_labels.size()) {
            label += quantum_number_l_labels.at(static_cast<size_t>(quantum_number_l_exp));
        } else {
            label += fmt::format("{:.0f}", quantum_number_l_exp);
        }
        label += "_";
    } else {
        label += fmt::format("S={:.1f},nu={:.1f},L={:.1f},F=", quantum_number_s_exp,
                             quantum_number_nu_exp, quantum_number_l_exp);
    }

    if (this->quantum_number_f == std::rintf(this->quantum_number_f)) {
        label += fmt::format("{:.0f}", this->quantum_number_f);
    } else if (2 * this->quantum_number_f == std::rintf(2 * this->quantum_number_f)) {
        label += fmt::format("{:.0f}/2", 2 * this->quantum_number_f);
    } else {
        std::abort(); // can't happen because the quantum number f is validated to be an integer
                      // or half-integer
    }

    if (this->quantum_number_m == std::rintf(this->quantum_number_m)) {
        label += fmt::format(",{:.0f}", this->quantum_number_m);
    } else if (2 * this->quantum_number_m == std::rintf(2 * this->quantum_number_m)) {
        label += fmt::format(",{:.0f}/2", 2 * this->quantum_number_m);
    } else {
        std::abort(); // can't happen because the quantum number m is validated to be an integer
                      // or half-integer
    }

    return label;
}

template <typename Real>
std::shared_ptr<KetAtom<Real>>
KetAtom<Real>::get_ket_for_different_quantum_number_m(Real new_quantum_number_m) const {
    auto ket = *this;
    ket.quantum_number_m = new_quantum_number_m;
    return std::make_shared<KetAtom<Real>>(ket);
}

template <typename Real>
const std::string &KetAtom<Real>::get_species() const {
    return species;
}

template <typename Real>
int KetAtom<Real>::get_quantum_number_n() const {
    return quantum_number_n;
}

template <typename Real>
Real KetAtom<Real>::get_quantum_number_nu() const {
    return quantum_number_nu_exp;
}

template <typename Real>
Real KetAtom<Real>::get_quantum_number_l() const {
    return quantum_number_l_exp;
}

template <typename Real>
Real KetAtom<Real>::get_quantum_number_s() const {
    return quantum_number_s_exp;
}

template <typename Real>
Real KetAtom<Real>::get_quantum_number_j() const {
    return quantum_number_j_exp;
}

template <typename Real>
Real KetAtom<Real>::get_quantum_number_nu_std() const {
    return quantum_number_nu_std;
}

template <typename Real>
Real KetAtom<Real>::get_quantum_number_l_std() const {
    return quantum_number_l_std;
}

template <typename Real>
Real KetAtom<Real>::get_quantum_number_s_std() const {
    return quantum_number_s_std;
}

template <typename Real>
Real KetAtom<Real>::get_quantum_number_j_std() const {
    return quantum_number_j_std;
}

template <typename Real>
bool KetAtom<Real>::operator==(const KetAtom<Real> &other) const {
    return Ket<Real>::operator==(other) && species == other.species &&
        quantum_number_n == other.quantum_number_n &&
        quantum_number_nu_exp == other.quantum_number_nu_exp &&
        quantum_number_nu_std == other.quantum_number_nu_std &&
        quantum_number_l_exp == other.quantum_number_l_exp &&
        quantum_number_l_std == other.quantum_number_l_std &&
        quantum_number_s_exp == other.quantum_number_s_exp &&
        quantum_number_s_std == other.quantum_number_s_std &&
        quantum_number_j_exp == other.quantum_number_j_exp &&
        quantum_number_j_std == other.quantum_number_j_std;
}

template <typename Real>
bool KetAtom<Real>::operator!=(const KetAtom<Real> &other) const {
    return !(*this == other);
}

template <typename Real>
size_t KetAtom<Real>::hash::operator()(const KetAtom<Real> &k) const {
    size_t seed = typename Ket<Real>::hash()(k);
    utils::hash_combine(seed, k.species);
    utils::hash_combine(seed, k.quantum_number_n);
    utils::hash_combine(seed, k.quantum_number_nu_exp);
    utils::hash_combine(seed, k.quantum_number_nu_std);
    utils::hash_combine(seed, k.quantum_number_l_exp);
    utils::hash_combine(seed, k.quantum_number_l_std);
    utils::hash_combine(seed, k.quantum_number_s_exp);
    utils::hash_combine(seed, k.quantum_number_s_std);
    utils::hash_combine(seed, k.quantum_number_j_exp);
    utils::hash_combine(seed, k.quantum_number_j_std);
    return seed;
}

// Explicit instantiations
template class KetAtom<float>;
template class KetAtom<double>;
} // namespace pairinteraction
