#pragma once

#include "pairinteraction/basis/Basis.hpp"
#include "pairinteraction/utils/Range.hpp"
#include "pairinteraction/utils/eigen_assertion.hpp"
#include "pairinteraction/utils/eigen_compat.hpp"
#include "pairinteraction/utils/hash.hpp"
#include "pairinteraction/utils/traits.hpp"

#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <complex>
#include <memory>
#include <unordered_map>
#include <vector>

namespace pairinteraction {
template <typename Scalar>
class BasisPairCreator;

template <typename Real>
class KetPair;

template <typename Scalar>
class BasisPair;

template <typename Scalar>
class BasisAtom;

template <typename Real>
class KetAtom;

template <typename Scalar>
struct traits::CrtpTraits<BasisPair<Scalar>> {
    using scalar_t = Scalar;
    using real_t = typename traits::NumTraits<Scalar>::real_t;
    using ket_t = KetPair<Scalar>;
    using ketvec_t = std::vector<std::shared_ptr<const ket_t>>;
};

template <typename Scalar>
class BasisPair : public Basis<BasisPair<Scalar>> {
    static_assert(traits::NumTraits<Scalar>::from_floating_point_v);

    friend class BasisPairCreator<Scalar>;
    struct Private {};

public:
    using Type = BasisPair<Scalar>;
    using real_t = typename traits::CrtpTraits<Type>::real_t;
    using ketvec_t = typename traits::CrtpTraits<Type>::ketvec_t;
    using range_t = Range<size_t>;
    using map_size_t = std::unordered_map<size_t, size_t>;
    using map_range_t = std::unordered_map<size_t, range_t>;
    using map_indices_t =
        std::unordered_map<std::vector<size_t>, size_t, utils::hash<std::vector<size_t>>>;

    BasisPair(Private /*unused*/, ketvec_t &&kets, map_range_t &&map_range_of_state_index2,
              map_indices_t &&state_indices_to_ket_index,
              std::shared_ptr<const BasisAtom<Scalar>> basis1,
              std::shared_ptr<const BasisAtom<Scalar>> basis2);
    const range_t &get_index_range(size_t state_index1) const;
    std::shared_ptr<const BasisAtom<Scalar>> get_basis1() const;
    std::shared_ptr<const BasisAtom<Scalar>> get_basis2() const;
    int get_ket_index_from_tuple(size_t state_index1, size_t state_index2) const;

    Eigen::VectorX<Scalar> get_amplitudes(std::shared_ptr<const KetAtom<real_t>> ket1,
                                          std::shared_ptr<const KetAtom<real_t>> ket2) const;
    Eigen::SparseMatrix<Scalar, Eigen::RowMajor>
    get_amplitudes(std::shared_ptr<const BasisAtom<Scalar>> other1,
                   std::shared_ptr<const BasisAtom<Scalar>> other2) const;
    Eigen::VectorX<real_t> get_overlaps(std::shared_ptr<const KetAtom<real_t>> ket1,
                                        std::shared_ptr<const KetAtom<real_t>> ket2) const;
    Eigen::SparseMatrix<real_t, Eigen::RowMajor>
    get_overlaps(std::shared_ptr<const BasisAtom<Scalar>> other1,
                 std::shared_ptr<const BasisAtom<Scalar>> other2) const;

private:
    map_range_t map_range_of_state_index2;
    map_indices_t state_indices_to_ket_index;
    std::shared_ptr<const BasisAtom<Scalar>> basis1;
    std::shared_ptr<const BasisAtom<Scalar>> basis2;
};

extern template class BasisPair<float>;
extern template class BasisPair<double>;
extern template class BasisPair<std::complex<float>>;
extern template class BasisPair<std::complex<double>>;
} // namespace pairinteraction
