// SPDX-FileCopyrightText: 2024 PairInteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include "pairinteraction/ket/Ket.hpp"
#include "pairinteraction/utils/traits.hpp"

#include <complex>
#include <initializer_list>
#include <memory>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <vector>

namespace pairinteraction {
template <typename Scalar>
class BasisPairCreator;

template <typename Scalar>
class BasisAtom;

/**
 * @class KetPair
 *
 * @brief Class for representing pair kets.
 *
 * Only quantum numbers that are well-defined for the pair state are contained in the map of
 * quantum numbers. The standard deviation of the contained quantum numbers is assumed to be zero.
 */
template <typename Scalar>
class KetPair : public Ket {
    static_assert(traits::NumTraits<Scalar>::from_floating_point_v);

    using real_t = typename traits::NumTraits<Scalar>::real_t;

    friend class BasisPairCreator<Scalar>;
    struct Private {};

public:
    KetPair(Private /*unused*/, std::initializer_list<size_t> atomic_indices,
            std::initializer_list<std::shared_ptr<const BasisAtom<Scalar>>> atomic_bases,
            real_t energy, std::unordered_map<std::string, double> quantum_numbers);

    bool has_quantum_number(const std::string &name) const;
    double get_quantum_number(const std::string &name) const;
    std::vector<std::shared_ptr<const BasisAtom<Scalar>>> get_atomic_states() const;

    bool operator==(const KetPair<Scalar> &other) const;
    bool operator!=(const KetPair<Scalar> &other) const;

    struct hash {
        std::size_t operator()(const KetPair<Scalar> &k) const;
    };

private:
    std::unordered_map<std::string, double> quantum_numbers;
    std::vector<size_t> atomic_indices;
    std::vector<std::shared_ptr<const BasisAtom<Scalar>>> atomic_bases;
};

extern template class KetPair<double>;
extern template class KetPair<std::complex<double>>;
} // namespace pairinteraction
