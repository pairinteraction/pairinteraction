#pragma once

#include "enums/SortBy.hpp"
#include "enums/TransformBy.hpp"
#include "utils/traits.hpp"

#include <Eigen/SparseCore>
#include <array>
#include <complex>
#include <vector>

template <typename Scalar>
class TransformableSortable {
public:
    using real_t = typename traits::NumTraits<Scalar>::real_t;

    virtual size_t get_number_of_states() const = 0;
    virtual size_t get_number_of_kets() const = 0;

    Eigen::SparseMatrix<Scalar> get_rotator(real_t alpha, real_t beta, real_t gamma) const;
    Eigen::SparseMatrix<Scalar> get_rotator(std::array<real_t, 3> to_z_axis,
                                            std::array<real_t, 3> to_y_axis) const;
    std::vector<int> get_sorter(SortBy label) const;
    std::vector<int> get_blocks(SortBy label) const;

    void transform(const Eigen::SparseMatrix<Scalar> &transformator);
    void rotate(real_t alpha, real_t beta, real_t gamma);
    void rotate(std::array<real_t, 3> to_z_axis, std::array<real_t, 3> to_y_axis);
    void sort(const std::vector<int> &sorter);
    void sort(SortBy label);

protected:
    virtual Eigen::SparseMatrix<Scalar> impl_get_rotator(real_t alpha, real_t beta,
                                                         real_t gamma) const = 0;
    virtual std::vector<int> impl_get_sorter(SortBy label) const = 0;
    virtual std::vector<int> impl_get_blocks(SortBy label) const = 0;

    virtual void impl_transform(const Eigen::SparseMatrix<Scalar> &transformator) = 0;
    virtual void impl_sort(const std::vector<int> &sorter) = 0;

    SortBy sorting{SortBy::KET};
    TransformBy transformation{TransformBy::IDENTITY};
};

extern template class TransformableSortable<float>;
extern template class TransformableSortable<double>;
extern template class TransformableSortable<std::complex<float>>;
extern template class TransformableSortable<std::complex<double>>;
