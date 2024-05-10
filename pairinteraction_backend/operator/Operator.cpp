#include "operator/Operator.hpp"
#include "basis/BasisAtom.hpp"
#include "enums/TransformationType.hpp"

#include <Eigen/SparseCore>
#include <memory>

template <typename Derived>
Operator<Derived>::Operator(std::shared_ptr<const basis_t> basis) : basis(basis) {}

template <typename Derived>
const Derived &Operator<Derived>::derived() const {
    return static_cast<const Derived &>(*this);
}

template <typename Derived>
std::shared_ptr<const typename Operator<Derived>::basis_t> Operator<Derived>::get_basis() const {
    return basis;
}

template <typename Derived>
const Eigen::SparseMatrix<typename Operator<Derived>::scalar_t, Eigen::RowMajor> &
Operator<Derived>::get_matrix() const {
    return matrix;
}

template <typename Derived>
const Transformation<typename Operator<Derived>::scalar_t> &
Operator<Derived>::get_transformation() const {
    return basis->get_transformation();
}

template <typename Derived>
Transformation<typename Operator<Derived>::scalar_t>
Operator<Derived>::get_rotator(real_t alpha, real_t beta, real_t gamma) const {
    return basis->get_rotator(alpha, beta, gamma);
}

template <typename Derived>
Sorting Operator<Derived>::get_sorter(TransformationType label) const {
    // Check that the label is a valid sorting label
    if (!utils::is_sorting(label)) {
        throw std::invalid_argument("The label is not a valid sorting label.");
    }

    // Get the sorter
    Sorting transformation =
        basis->get_sorter_without_checks(label & ~TransformationType::SORT_BY_ENERGY);

    if ((label & TransformationType::SORT_BY_ENERGY) == TransformationType::SORT_BY_ENERGY) {
        std::vector<real_t> energies_of_states;
        energies_of_states.reserve(matrix.rows());
        for (int i = 0; i < matrix.rows(); ++i) {
            if constexpr (traits::is_complex_v<scalar_t>) {
                energies_of_states.push_back(matrix.coeff(i, i).real());
            } else {
                energies_of_states.push_back(matrix.coeff(i, i));
            }
        }

        std::stable_sort(
            transformation.matrix.indices().data(),
            transformation.matrix.indices().data() + transformation.matrix.indices().size(),
            [&](int i, int j) { return energies_of_states[i] < energies_of_states[j]; });

        transformation.transformation_type.push_back(TransformationType::SORT_BY_ENERGY);
    }

    // Check if the full label has been used for sorting
    if (!utils::is_comprised_by_label(label, transformation.transformation_type)) {
        throw std::invalid_argument("The states could not be sorted by the requested label.");
    }

    return transformation;
}

template <typename Derived>
Blocks Operator<Derived>::get_blocks(TransformationType label) const {
    // Check that the label is a valid sorting label
    if (!utils::is_sorting(label)) {
        throw std::invalid_argument("The label is not a valid sorting label.");
    }

    // Check if the states are sorted by the requested label
    if (!utils::is_sorted_by_label(label, get_transformation().transformation_type)) {
        throw std::invalid_argument("The states are not sorted by the requested label.");
    }

    // Get the blocks
    Blocks blocks = basis->get_blocks_without_checks(label & ~TransformationType::SORT_BY_ENERGY);

    if ((label & TransformationType::SORT_BY_ENERGY) == TransformationType::SORT_BY_ENERGY) {
        std::vector<int> completed;
        size_t block_idx = 0;
        scalar_t last_diagonal = matrix.coeff(0, 0);

        for (int i = 0; i < matrix.rows(); ++i) {
            if (matrix.coeff(i, i) != last_diagonal) {
                completed.push_back(i);
            } else if (block_idx < blocks.start.size() && i == blocks.start[block_idx]) {
                completed.push_back(i);
                ++block_idx;
            }
        }

        blocks.start = completed;
        blocks.transformation_type.push_back(TransformationType::SORT_BY_ENERGY);
    }

    // Check if the full label has been used for getting the blocks
    if (!utils::is_comprised_by_label(label, blocks.transformation_type)) {
        throw std::invalid_argument("The blocks could not be obtained by the requested label.");
    }

    return blocks;
}

template <typename Derived>
Derived Operator<Derived>::transform(
    const Transformation<typename Operator<Derived>::scalar_t> &transformation) const {
    auto transformed = derived();
    transformed.matrix = transformation.matrix.adjoint() * matrix * transformation.matrix;
    transformed.basis = basis->transform(transformation);
    return transformed;
}

template <typename Derived>
Derived Operator<Derived>::transform(const Sorting &transformation) const {
    auto transformed = derived();
    transformed.matrix = matrix.twistedBy(transformation.matrix.inverse());
    transformed.basis = basis->transform(transformation);
    return transformed;
}

// Overloaded operators
template <typename Derived>
Derived operator*(const typename Operator<Derived>::scalar_t &lhs, const Operator<Derived> &rhs) {
    Derived result = rhs.derived();
    result.matrix *= lhs;
    return result;
}

template <typename Derived>
Derived operator*(const Operator<Derived> &lhs, const typename Operator<Derived>::scalar_t &rhs) {
    Derived result = lhs.derived();
    result.matrix *= rhs;
    return result;
}

template <typename Derived>
Derived operator/(const Operator<Derived> &lhs, const typename Operator<Derived>::scalar_t &rhs) {
    Derived result = lhs.derived();
    result.matrix /= rhs;
    return result;
}

template <typename Derived>
Derived operator+(const Operator<Derived> &lhs, const Operator<Derived> &rhs) {
    if (lhs.basis != rhs.basis) {
        throw std::invalid_argument("The basis of the operators is not the same.");
    }
    Derived result = lhs.derived();
    result.matrix += rhs.matrix;
    return result;
}

template <typename Derived>
Derived operator-(const Operator<Derived> &lhs, const Operator<Derived> &rhs) {
    if (lhs.basis != rhs.basis) {
        throw std::invalid_argument("The basis of the operators is not the same.");
    }
    Derived result = lhs.derived();
    result.matrix -= rhs.matrix;
    return result;
}

// Explicit instantiations
#include "operator/OperatorAtom.hpp"

template class Operator<OperatorAtom<float>>;
template class Operator<OperatorAtom<double>>;
template class Operator<OperatorAtom<std::complex<float>>>;
template class Operator<OperatorAtom<std::complex<double>>>;

template OperatorAtom<float> operator*(const float &lhs, const Operator<OperatorAtom<float>> &rhs);
template OperatorAtom<double> operator*(const double &lhs,
                                        const Operator<OperatorAtom<double>> &rhs);
template OperatorAtom<std::complex<float>>
operator*(const std::complex<float> &lhs, const Operator<OperatorAtom<std::complex<float>>> &rhs);
template OperatorAtom<std::complex<double>>
operator*(const std::complex<double> &lhs, const Operator<OperatorAtom<std::complex<double>>> &rhs);

template OperatorAtom<float> operator*(const Operator<OperatorAtom<float>> &lhs, const float &rhs);
template OperatorAtom<double> operator*(const Operator<OperatorAtom<double>> &lhs,
                                        const double &rhs);
template OperatorAtom<std::complex<float>>
operator*(const Operator<OperatorAtom<std::complex<float>>> &lhs, const std::complex<float> &rhs);
template OperatorAtom<std::complex<double>>
operator*(const Operator<OperatorAtom<std::complex<double>>> &lhs, const std::complex<double> &rhs);

template OperatorAtom<float> operator/(const Operator<OperatorAtom<float>> &lhs, const float &rhs);
template OperatorAtom<double> operator/(const Operator<OperatorAtom<double>> &lhs,
                                        const double &rhs);
template OperatorAtom<std::complex<float>>
operator/(const Operator<OperatorAtom<std::complex<float>>> &lhs, const std::complex<float> &rhs);
template OperatorAtom<std::complex<double>>
operator/(const Operator<OperatorAtom<std::complex<double>>> &lhs, const std::complex<double> &rhs);

template OperatorAtom<float> operator+(const Operator<OperatorAtom<float>> &lhs,
                                       const Operator<OperatorAtom<float>> &rhs);
template OperatorAtom<double> operator+(const Operator<OperatorAtom<double>> &lhs,
                                        const Operator<OperatorAtom<double>> &rhs);
template OperatorAtom<std::complex<float>>
operator+(const Operator<OperatorAtom<std::complex<float>>> &lhs,
          const Operator<OperatorAtom<std::complex<float>>> &rhs);
template OperatorAtom<std::complex<double>>
operator+(const Operator<OperatorAtom<std::complex<double>>> &lhs,
          const Operator<OperatorAtom<std::complex<double>>> &rhs);

template OperatorAtom<float> operator-(const Operator<OperatorAtom<float>> &lhs,
                                       const Operator<OperatorAtom<float>> &rhs);
template OperatorAtom<double> operator-(const Operator<OperatorAtom<double>> &lhs,
                                        const Operator<OperatorAtom<double>> &rhs);
template OperatorAtom<std::complex<float>>
operator-(const Operator<OperatorAtom<std::complex<float>>> &lhs,
          const Operator<OperatorAtom<std::complex<float>>> &rhs);
template OperatorAtom<std::complex<double>>
operator-(const Operator<OperatorAtom<std::complex<double>>> &lhs,
          const Operator<OperatorAtom<std::complex<double>>> &rhs);
