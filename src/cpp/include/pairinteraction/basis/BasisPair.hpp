// SPDX-FileCopyrightText: 2024 PairInteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include "pairinteraction/basis/Basis.hpp"
#include "pairinteraction/enums/SorterType.hpp"
#include "pairinteraction/utils/Range.hpp"
#include "pairinteraction/utils/eigen_assertion.hpp"
#include "pairinteraction/utils/eigen_compat.hpp"
#include "pairinteraction/utils/hash.hpp"
#include "pairinteraction/utils/traits.hpp"

#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <complex>
#include <map>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

namespace pairinteraction {
enum class OperatorType;

template <typename Scalar>
class BasisPairCreator;

template <typename Scalar>
class KetPair;

template <typename Scalar>
class BasisPair;

template <typename Scalar>
class BasisAtom;

template <typename Scalar>
struct traits::CrtpTraits<BasisPair<Scalar>> {
    using scalar_t = Scalar;
    using real_t = typename traits::NumTraits<Scalar>::real_t;
    using ket_t = KetPair<Scalar>;
    using ketvec_t = std::vector<std::shared_ptr<const ket_t>>;
};

template <typename Scalar>
class BasisPair : public Basis<BasisPair<Scalar>>,
                  public std::enable_shared_from_this<BasisPair<Scalar>> {
    static_assert(traits::NumTraits<Scalar>::from_floating_point_v);

    friend class BasisPairCreator<Scalar>;
    struct Private {};

public:
    using Type = BasisPair<Scalar>;
    using real_t = typename traits::CrtpTraits<Type>::real_t;
    using quantum_numbers_t = typename Basis<Type>::quantum_numbers_t;
    using energies_t = typename Basis<Type>::energies_t;
    using ket_t = typename traits::CrtpTraits<Type>::ket_t;
    using ketvec_t = typename traits::CrtpTraits<Type>::ketvec_t;
    using range_t = Range<size_t>;
    using map_size_t = std::unordered_map<size_t, size_t>;
    using map_range_t = std::unordered_map<size_t, range_t>;
    using map_indices_t =
        std::unordered_map<std::vector<size_t>, size_t, utils::hash<std::vector<size_t>>>;

    static const std::map<SorterType, std::string> sorter_type_to_quantum_number_name;

    BasisPair(Private /*unused*/, std::shared_ptr<const ketvec_t> kets,
              std::shared_ptr<const map_range_t> state_index1_to_state_index_range2,
              std::shared_ptr<const map_indices_t> state_indices_to_ket_index,
              std::shared_ptr<const BasisAtom<Scalar>> basis1,
              std::shared_ptr<const BasisAtom<Scalar>> basis2,
              Eigen::SparseMatrix<Scalar, Eigen::RowMajor> &&coefficients,
              energies_t &&energy_of_states, quantum_numbers_t &&quantum_numbers_of_states);
    const range_t &get_index_range(size_t state_index1) const;
    std::shared_ptr<const BasisAtom<Scalar>> get_basis1() const;
    std::shared_ptr<const BasisAtom<Scalar>> get_basis2() const;
    int get_ket_index_from_tuple(size_t state_index1, size_t state_index2) const;
    std::shared_ptr<const Type> merge(const std::shared_ptr<const Type> &other) const;
    // Create a copy of the basis that has the same kets but new states. Quantum numbers that are
    // not given are not well-defined for the states, and so is the energy if it is not given.
    // Given energies must be finite and given quantum numbers must be integers or half-integers;
    // the maximum value is accepted for neither.
    std::shared_ptr<const Type>
    get_modified_copy(Eigen::SparseMatrix<Scalar, Eigen::RowMajor> coefficients,
                      energies_t energy_of_states,
                      quantum_numbers_t quantum_numbers_of_states) const;

    Eigen::SparseMatrix<Scalar, Eigen::RowMajor>
    get_matrix_elements(const std::shared_ptr<const Type> &final_state, OperatorType type1,
                        OperatorType type2, int q1, int q2) const;

private:
    // The following maps refer to the kets and are thus shared like them
    // state index 1 -> (min state index 2, max state index 2)
    std::shared_ptr<const map_range_t> state_index1_to_state_index_range2;
    // (state atom index 1, state atom index 2) -> ket pair index
    std::shared_ptr<const map_indices_t> state_indices_to_ket_index;
    std::shared_ptr<const BasisAtom<Scalar>> basis1;
    std::shared_ptr<const BasisAtom<Scalar>> basis2;
};

extern template class BasisPair<double>;
extern template class BasisPair<std::complex<double>>;
} // namespace pairinteraction
