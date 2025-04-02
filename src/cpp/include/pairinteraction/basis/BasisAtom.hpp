// SPDX-FileCopyrightText: 2024 Pairinteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include "pairinteraction/basis/Basis.hpp"
#include "pairinteraction/utils/traits.hpp"

#include <complex>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

namespace pairinteraction {
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

    BasisAtom(Private /*unused*/, ketvec_t &&kets, std::string &&id_of_kets, Database &database);
    Database &get_database() const;
    const std::string &get_species() const;
    const std::string &get_id_of_kets() const;

    int get_ket_index_from_id(size_t ket_id) const;

    Eigen::VectorX<Scalar> get_matrix_elements(std::shared_ptr<const ket_t> ket, OperatorType type,
                                               int q = 0) const override;
    Eigen::SparseMatrix<Scalar, Eigen::RowMajor>
    get_matrix_elements(std::shared_ptr<const Type> other, OperatorType type,
                        int q = 0) const override;

private:
    std::string id_of_kets;
    Database &database;
    std::string species;
    std::unordered_map<size_t, size_t> ket_id_to_ket_index;
};

extern template class BasisAtom<double>;
extern template class BasisAtom<std::complex<double>>;
} // namespace pairinteraction
