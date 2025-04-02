// SPDX-FileCopyrightText: 2024 Pairinteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include "pairinteraction/operator/Operator.hpp"
#include "pairinteraction/utils/traits.hpp"

#include <memory>
#include <vector>

namespace pairinteraction {
class Database;

enum class OperatorType;

template <typename Scalar>
class BasisAtom;

class KetAtom;

template <typename T>
class OperatorAtom;

template <typename Scalar>
struct traits::CrtpTraits<OperatorAtom<Scalar>> {
    using scalar_t = Scalar;
    using real_t = typename traits::NumTraits<Scalar>::real_t;
    using ket_t = KetAtom;
    using ketvec_t = std::vector<std::shared_ptr<const ket_t>>;
    using basis_t = BasisAtom<scalar_t>;
};

template <typename Scalar>
class OperatorAtom : public Operator<OperatorAtom<Scalar>> {
public:
    static_assert(traits::NumTraits<Scalar>::from_floating_point_v);

    using Type = OperatorAtom<Scalar>;
    using basis_t = typename traits::CrtpTraits<Type>::basis_t;

    OperatorAtom(std::shared_ptr<const basis_t> basis);
    OperatorAtom(std::shared_ptr<const basis_t> basis, OperatorType type, int q = 0);

private:
    friend class Database;

    OperatorAtom(std::shared_ptr<const basis_t> basis,
                 Eigen::SparseMatrix<Scalar, Eigen::RowMajor> &&matrix);
};

extern template class OperatorAtom<double>;
extern template class OperatorAtom<std::complex<double>>;
} // namespace pairinteraction
