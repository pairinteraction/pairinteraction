#pragma once

#include "pairinteraction/interfaces/TransformationBuilderInterface.hpp"
#include "pairinteraction/utils/eigen_assertion.hpp"
#include "pairinteraction/utils/eigen_compat.hpp"
#include "pairinteraction/utils/traits.hpp"

#include <Eigen/SparseCore>
#include <memory>
#include <set>
#include <vector>

namespace pairinteraction {
enum class TransformationType : unsigned char;

template <typename Scalar>
class DiagonalizerInterface;

template <typename Sortable>
class Range;

template <typename Derived>
class System
    : public TransformationBuilderInterface<typename traits::CrtpTraits<Derived>::scalar_t> {
public:
    using scalar_t = typename traits::CrtpTraits<Derived>::scalar_t;
    using real_t = typename traits::CrtpTraits<Derived>::real_t;
    using ketvec_t = typename traits::CrtpTraits<Derived>::ketvec_t;
    using basis_t = typename traits::CrtpTraits<Derived>::basis_t;
    using operator_t = typename traits::CrtpTraits<Derived>::operator_t;

    System(std::shared_ptr<const basis_t> basis);
    System(const System &other);
    System(System &&other) noexcept;
    System<Derived> &operator=(const System &other);
    System<Derived> &operator=(System &&other) noexcept;
    virtual ~System();

    std::shared_ptr<const basis_t> get_basis() const;

    const Eigen::SparseMatrix<scalar_t, Eigen::RowMajor> &get_matrix() const;

    const Transformation<scalar_t> &get_transformation() const override;
    Transformation<scalar_t> get_rotator(real_t alpha, real_t beta, real_t gamma) const override;
    Sorting get_sorter(const std::vector<TransformationType> &labels) const override;
    std::vector<IndicesOfBlock>
    get_indices_of_blocks(const std::vector<TransformationType> &labels) const override;

    void get_indices_of_blocks_without_checks(const std::set<TransformationType> &unique_labels,
                                              IndicesOfBlocksCreator &blocks) const;

    Derived transformed(const Transformation<scalar_t> &transformation) const;
    Derived transformed(const Sorting &transformation) const;

    System<Derived> &diagonalize(const DiagonalizerInterface<scalar_t> &diagonalizer,
                                 int precision = 12, const Range<real_t> &eigenvalue_range = {});
    bool is_diagonal() const;
    const Eigen::SparseMatrix<scalar_t, Eigen::RowMajor> &get_eigenstates() const;
    Eigen::VectorX<real_t> get_eigenvalues() const;

protected:
    mutable std::unique_ptr<operator_t> hamiltonian;
    mutable bool hamiltonian_requires_construction{true};
    mutable bool hamiltonian_is_diagonal{false};
    mutable std::vector<TransformationType> blockdiagonalizing_labels;

    virtual void construct_hamiltonian() const = 0;

private:
    const Derived &derived() const;
};
} // namespace pairinteraction
