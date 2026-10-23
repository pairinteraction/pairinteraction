// SPDX-FileCopyrightText: 2024 PairInteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include "pairinteraction/basis/Basis.hpp"
#include "pairinteraction/enums/SorterType.hpp"
#include "pairinteraction/utils/traits.hpp"

#include <complex>
#include <map>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

namespace pairinteraction {
enum class OperatorType;

class Database;

class KetAtom;

template <typename Scalar>
class BasisAtom;

template <typename Scalar>
struct traits::CrtpTraits<BasisAtom<Scalar>> {
    using scalar_t = Scalar;
    using real_t = typename traits::NumTraits<Scalar>::real_t;
    using ket_t = KetAtom;
    using ketvec_t = std::vector<std::shared_ptr<const ket_t>>;
};

/**
 * @class BasisAtom
 *
 * @brief Class for creating a basis of atomic kets.
 *
 * @tparam Scalar Complex number type.
 */
template <typename Scalar>
class BasisAtom : public Basis<BasisAtom<Scalar>>,
                  public std::enable_shared_from_this<BasisAtom<Scalar>> {
    static_assert(traits::NumTraits<Scalar>::from_floating_point_v);

    friend class Database;
    struct Private {};

public:
    using Type = BasisAtom<Scalar>;
    using ket_t = typename traits::CrtpTraits<Type>::ket_t;
    using ketvec_t = typename traits::CrtpTraits<Type>::ketvec_t;
    using real_t = typename traits::CrtpTraits<Type>::real_t;
    using quantum_numbers_t = typename Basis<Type>::quantum_numbers_t;
    using energies_t = typename Basis<Type>::energies_t;
    using map_size_t = std::unordered_map<size_t, size_t>;

    static const std::map<SorterType, std::string> sorter_type_to_quantum_number_name;

    BasisAtom(Private /*unused*/, ketvec_t &&kets, map_size_t &&ket_id_to_ket_index,
              std::string &&canonical_basis_id, Database &database,
              Eigen::SparseMatrix<Scalar, Eigen::RowMajor> &&coefficients,
              energies_t &&energy_of_states, quantum_numbers_t &&quantum_numbers_of_states);
    Database &get_database() const;
    const std::string &get_species() const;
    const std::string &get_canonical_basis_id() const;
    std::shared_ptr<const Type> merge(std::shared_ptr<const Type> other) const override;
    // Create a copy of the basis that has the same kets but new states. Quantum numbers that are
    // not given are not well-defined, and so are quantum numbers and energies whose value is the
    // maximum.
    std::shared_ptr<const Type>
    get_modified_copy(Eigen::SparseMatrix<Scalar, Eigen::RowMajor> coefficients,
                      energies_t energy_of_states,
                      quantum_numbers_t quantum_numbers_of_states) const;

    int get_ket_index_from_id(size_t ket_id) const;

    Eigen::SparseMatrix<Scalar, Eigen::RowMajor>
    get_matrix_elements(const std::shared_ptr<const Type> &other, OperatorType type, int q) const;

private:
    std::string canonical_basis_id;
    Database &database;
    // ket id in database -> ket index
    map_size_t ket_id_to_ket_index;
};

extern template class BasisAtom<double>;
extern template class BasisAtom<std::complex<double>>;
} // namespace pairinteraction
