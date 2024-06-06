#pragma once

#include <memory>

#include "interfaces/TransformationBuilderInterface.hpp"
#include "utils/traits.hpp"

enum class TransformationType : unsigned char;

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
    virtual ~System();

    std::shared_ptr<const basis_t> get_basis() const;

    const Eigen::SparseMatrix<scalar_t, Eigen::RowMajor> &get_matrix() const;

    const Transformation<scalar_t> &get_transformation() const override;
    Transformation<scalar_t> get_rotator(real_t alpha, real_t beta, real_t gamma) const override;
    Sorting get_sorter(TransformationType label) const override;
    Blocks get_blocks(TransformationType label) const override;

    Derived transformed(const Transformation<scalar_t> &transformation) const;
    Derived transformed(const Sorting &transformation) const;

protected:
    mutable std::unique_ptr<operator_t> hamiltonian;
    mutable bool hamiltonian_requires_construction{true};

    virtual void construct_hamiltonian() const = 0;

private:
    const Derived &derived() const;
};
