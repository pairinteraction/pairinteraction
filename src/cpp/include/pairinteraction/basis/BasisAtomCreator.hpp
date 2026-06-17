// SPDX-FileCopyrightText: 2024 PairInteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include "pairinteraction/utils/Range.hpp"
#include "pairinteraction/utils/traits.hpp"

#include <complex>
#include <memory>
#include <optional>
#include <string>
#include <unordered_map>
#include <vector>

namespace pairinteraction {
template <typename Scalar>
class BasisAtom;

class KetAtom;

class Database;

/**
 * @class BasisAtomCreator
 *
 * @brief Builder class for creating BasisAtom objects.
 *
 * @tparam Scalar Complex number type.
 */
template <typename Scalar>
class BasisAtomCreator {
    static_assert(traits::NumTraits<Scalar>::from_floating_point_v);

public:
    using real_t = typename traits::NumTraits<Scalar>::real_t;
    using ket_t = KetAtom;
    BasisAtomCreator() = default;
    BasisAtomCreator<Scalar> &set_species(const std::string &value);
    BasisAtomCreator<Scalar> &restrict_energy(real_t min, real_t max);
    // Set the quantum number range with the given logical name (e.g. "f", "m", "n", "l", ...).
    // The parity is restricted via the name "parity" with a value of +1 (even) or -1 (odd).
    BasisAtomCreator<Scalar> &restrict_quantum_number(const std::string &name, real_t min,
                                                      real_t max);
    BasisAtomCreator<Scalar> &set_quantum_number_standard_deviation_factor(real_t value);
    BasisAtomCreator<Scalar> &add_ket(const std::shared_ptr<const ket_t> &ket);
    std::shared_ptr<const BasisAtom<Scalar>> create(Database &database) const;

private:
    std::optional<std::string> species;
    Range<real_t> range_energy;
    std::unordered_map<std::string, Range<real_t>> quantum_number_ranges;
    real_t quantum_number_standard_deviation_factor{2};
    std::vector<size_t> additional_ket_ids;
    std::optional<std::string> additional_ket_species;
};

extern template class BasisAtomCreator<double>;
extern template class BasisAtomCreator<std::complex<double>>;
} // namespace pairinteraction
