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
class BasisPair;

template <typename Scalar>
class BasisAtom;

template <typename Scalar>
class KetPair;

template <typename Scalar>
class BasisPairCreator {
    static_assert(traits::NumTraits<Scalar>::from_floating_point_v);

public:
    using real_t = typename traits::NumTraits<Scalar>::real_t;
    using basis_t = BasisPair<Scalar>;
    using ket_t = KetPair<Scalar>;
    using ketvec_t = std::vector<std::shared_ptr<const ket_t>>;

    BasisPairCreator(std::shared_ptr<const BasisAtom<Scalar>> basis1,
                     std::shared_ptr<const BasisAtom<Scalar>> basis2);
    BasisPairCreator<Scalar> &restrict_energy(real_t min, real_t max);
    // Set the quantum number range with the given name ("m", "parity_under_inversion", ...).
    // The parities are restricted with values of +1 (even) or -1 (odd).
    BasisPairCreator<Scalar> &restrict_quantum_number(const std::string &name, real_t min,
                                                      real_t max);
    BasisPairCreator<Scalar> &set_symmetrization_enabled(bool value);
    std::shared_ptr<const BasisPair<Scalar>> create() const;

private:
    std::shared_ptr<const BasisAtom<Scalar>> basis1;
    std::shared_ptr<const BasisAtom<Scalar>> basis2;
    Range<real_t> range_energy;
    std::unordered_map<std::string, Range<real_t>> quantum_number_ranges;
    bool symmetrization_enabled;
};

extern template class BasisPairCreator<double>;
extern template class BasisPairCreator<std::complex<double>>;
} // namespace pairinteraction
