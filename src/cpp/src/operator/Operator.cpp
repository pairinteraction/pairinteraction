#include "pairinteraction/operator/Operator.hpp"

#include "pairinteraction/basis/BasisAtom.hpp"
#include "pairinteraction/enums/TransformationType.hpp"
#include "pairinteraction/operator/OperatorAtom.hpp"
#include "pairinteraction/utils/eigen_assertion.hpp"

#include <Eigen/SparseCore>
#include <memory>

namespace pairinteraction {
template <typename Derived>
Operator<Derived>::Operator(std::shared_ptr<const basis_t> basis) : basis(std::move(basis)) {}

template <typename Derived>
const Derived &Operator<Derived>::derived() const {
    return static_cast<const Derived &>(*this);
}

template <typename Derived>
Derived &Operator<Derived>::derived_mutable() {
    return static_cast<Derived &>(*this);
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
Sorting Operator<Derived>::get_sorter(const std::vector<TransformationType> &labels) const {
    basis->perform_sorter_checks(labels);

    // Split labels into three parts (one before SORT_BY_ENERGY, one with SORT_BY_ENERGY, and one
    // after)
    auto it = std::find(labels.begin(), labels.end(), TransformationType::SORT_BY_ENERGY);
    std::vector<TransformationType> before_energy(labels.begin(), it);
    bool contains_energy = (it != labels.end());
    std::vector<TransformationType> after_energy(contains_energy ? it + 1 : labels.end(),
                                                 labels.end());

    // Initialize transformation
    Sorting transformation;

    // Apply sorting for labels before SORT_BY_ENERGY
    if (!before_energy.empty()) {
        basis->get_sorter_without_checks(before_energy, transformation);
    }

    // Apply SORT_BY_ENERGY if present
    if (contains_energy) {
        std::vector<real_t> energies_of_states;
        energies_of_states.reserve(matrix.rows());
        for (int i = 0; i < matrix.rows(); ++i) {
            energies_of_states.push_back(std::real(matrix.coeff(i, i)));
        }

        if (transformation.matrix.size() == 0) {
            // If no previous sorting, create identity permutation
            transformation.matrix.resize(matrix.rows());
            transformation.matrix.setIdentity();
        }

        std::stable_sort(
            transformation.matrix.indices().data(),
            transformation.matrix.indices().data() + transformation.matrix.indices().size(),
            [&](int i, int j) { return energies_of_states[i] < energies_of_states[j]; });

        transformation.transformation_type.push_back(TransformationType::SORT_BY_ENERGY);
    }

    // Apply sorting for labels after SORT_BY_ENERGY
    if (!after_energy.empty()) {
        basis->get_sorter_without_checks(after_energy, transformation);
    }

    // Check if all labels have been used for sorting
    if (!utils::are_same_labels(labels, transformation.transformation_type)) {
        throw std::invalid_argument("The states could not be sorted by all the requested labels.");
    }

    return transformation;
}

template <typename Derived>
Blocks Operator<Derived>::get_blocks(const std::vector<TransformationType> &labels) const {
    std::set<TransformationType> unique_labels(labels.begin(), labels.end());
    basis->perform_blocks_checks(unique_labels);

    // Split labels into two parts (one with SORT_BY_ENERGY and one without)
    auto it = unique_labels.find(TransformationType::SORT_BY_ENERGY);
    bool contains_energy = (it != unique_labels.end());
    if (contains_energy) {
        unique_labels.erase(it);
    }

    // Initialize blocks
    Blocks blocks;

    // Handle all labels except SORT_BY_ENERGY
    if (!unique_labels.empty()) {
        basis->get_blocks_without_checks(unique_labels, blocks);
    }

    // Handle SORT_BY_ENERGY if present
    if (contains_energy) {
        scalar_t last_energy = std::real(matrix.coeff(0, 0));

        std::vector<int> blocks_start_original;
        blocks_start_original.swap(blocks.start);
        blocks.start.reserve(matrix.rows());
        size_t idx = 0;

        for (int i = 0; i < matrix.rows(); ++i) {
            if (idx < blocks_start_original.size() && i == blocks_start_original[idx]) {
                blocks.start.push_back(i);
                ++idx;
            } else if (std::real(matrix.coeff(i, i)) != last_energy) {
                blocks.start.push_back(i);
                last_energy = std::real(matrix.coeff(i, i));
            }
        }

        blocks.start.shrink_to_fit();
    }

    return blocks;
}

template <typename Derived>
Derived Operator<Derived>::transformed(
    const Transformation<typename Operator<Derived>::scalar_t> &transformation) const {
    auto transformed = derived();
    transformed.matrix = transformation.matrix.adjoint() * matrix * transformation.matrix;
    transformed.basis = basis->transformed(transformation);
    return transformed;
}

template <typename Derived>
Derived Operator<Derived>::transformed(const Sorting &transformation) const {
    auto transformed = derived();
    transformed.matrix = matrix.twistedBy(transformation.matrix.inverse());
    transformed.basis = basis->transformed(transformation);
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
Derived &operator+=(Operator<Derived> &lhs, const Operator<Derived> &rhs) {
    if (lhs.basis != rhs.basis) {
        throw std::invalid_argument("The basis of the operators is not the same.");
    }
    lhs.matrix += rhs.matrix;
    return lhs.derived_mutable();
}

template <typename Derived>
Derived &operator-=(Operator<Derived> &lhs, const Operator<Derived> &rhs) {
    if (lhs.basis != rhs.basis) {
        throw std::invalid_argument("The basis of the operators is not the same.");
    }
    lhs.matrix -= rhs.matrix;
    return lhs.derived_mutable();
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

template OperatorAtom<float> &operator+=(Operator<OperatorAtom<float>> &lhs,
                                         const Operator<OperatorAtom<float>> &rhs);
template OperatorAtom<double> &operator+=(Operator<OperatorAtom<double>> &lhs,
                                          const Operator<OperatorAtom<double>> &rhs);
template OperatorAtom<std::complex<float>> &
operator+=(Operator<OperatorAtom<std::complex<float>>> &lhs,
           const Operator<OperatorAtom<std::complex<float>>> &rhs);
template OperatorAtom<std::complex<double>> &
operator+=(Operator<OperatorAtom<std::complex<double>>> &lhs,
           const Operator<OperatorAtom<std::complex<double>>> &rhs);

template OperatorAtom<float> &operator-=(Operator<OperatorAtom<float>> &lhs,
                                         const Operator<OperatorAtom<float>> &rhs);
template OperatorAtom<double> &operator-=(Operator<OperatorAtom<double>> &lhs,
                                          const Operator<OperatorAtom<double>> &rhs);
template OperatorAtom<std::complex<float>> &
operator-=(Operator<OperatorAtom<std::complex<float>>> &lhs,
           const Operator<OperatorAtom<std::complex<float>>> &rhs);
template OperatorAtom<std::complex<double>> &
operator-=(Operator<OperatorAtom<std::complex<double>>> &lhs,
           const Operator<OperatorAtom<std::complex<double>>> &rhs);

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
} // namespace pairinteraction
