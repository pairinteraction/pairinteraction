// SPDX-FileCopyrightText: 2024 Pairinteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include "pairinteraction/interfaces/TransformationBuilderInterface.hpp"
#include "pairinteraction/utils/eigen_assertion.hpp"
#include "pairinteraction/utils/eigen_compat.hpp"
#include "pairinteraction/utils/traits.hpp"

#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <memory>
#include <set>
#include <unordered_map>
#include <vector>

namespace pairinteraction {
enum class Parity : int;
enum class TransformationType : unsigned char;
enum class OperatorType;

/**
 * @class Basis
 *
 * @brief Base class for a basis
 *
 * This base class represents a basis. It comprises a list of ket states and a matrix of
 * coefficients. The rows of the coefficient matrix correspond to indices of ket states and
 * the columns to indices of basis vectors.
 * Using CRPT, it is a base class for specific basis implementations. Its
 * constructor is protected to indicate that derived classes should not allow direct instantiation.
 * Instead, a factory class should be provided that is a friend of the derived class and can create
 * instances of it.
 *
 * @tparam Derived Derived class.
 */

template <typename Derived>
class Basis
    : public TransformationBuilderInterface<typename traits::CrtpTraits<Derived>::scalar_t> {
public:
    using scalar_t = typename traits::CrtpTraits<Derived>::scalar_t;
    using real_t = typename traits::CrtpTraits<Derived>::real_t;
    using ket_t = typename traits::CrtpTraits<Derived>::ket_t;
    using ketvec_t = typename traits::CrtpTraits<Derived>::ketvec_t;

    Basis() = delete;
    virtual ~Basis() = default;

    bool has_quantum_number_f() const;
    bool has_quantum_number_m() const;
    bool has_parity() const;

    const ketvec_t &get_kets() const;
    size_t get_number_of_states() const;
    size_t get_number_of_kets() const;
    real_t get_quantum_number_f(size_t state_index) const;
    real_t get_quantum_number_m(size_t state_index) const;
    Parity get_parity(size_t state_index) const;
    std::shared_ptr<const Derived> get_state(size_t state_index) const;
    std::shared_ptr<const ket_t> get_ket(size_t ket_index) const;
    std::shared_ptr<const ket_t> get_corresponding_ket(size_t state_index) const;
    std::shared_ptr<const ket_t> get_corresponding_ket(std::shared_ptr<const Derived> state) const;
    size_t get_corresponding_ket_index(size_t state_index) const;
    size_t get_corresponding_ket_index(std::shared_ptr<const Derived> state) const;
    std::shared_ptr<const Derived> get_corresponding_state(size_t ket_index) const;
    std::shared_ptr<const Derived> get_corresponding_state(std::shared_ptr<const ket_t> ket) const;
    size_t get_corresponding_state_index(size_t ket_index) const;
    size_t get_corresponding_state_index(std::shared_ptr<const ket_t> ket) const;
    std::shared_ptr<const Derived> get_canonical_state_from_ket(size_t ket_index) const;
    std::shared_ptr<const Derived>
    get_canonical_state_from_ket(std::shared_ptr<const ket_t> ket) const;
    const Eigen::SparseMatrix<scalar_t, Eigen::RowMajor> &get_coefficients() const;
    Eigen::SparseMatrix<scalar_t, Eigen::RowMajor> &get_coefficients();
    Eigen::VectorX<scalar_t> get_amplitudes(std::shared_ptr<const ket_t> ket) const;
    Eigen::SparseMatrix<scalar_t, Eigen::RowMajor>
    get_amplitudes(std::shared_ptr<const Derived> other) const;
    Eigen::VectorX<real_t> get_overlaps(std::shared_ptr<const ket_t> ket) const;
    Eigen::SparseMatrix<real_t, Eigen::RowMajor>
    get_overlaps(std::shared_ptr<const Derived> other) const;
    virtual Eigen::VectorX<scalar_t> get_matrix_elements(std::shared_ptr<const ket_t> ket,
                                                         OperatorType type, int q = 0) const = 0;
    Eigen::SparseMatrix<scalar_t, Eigen::RowMajor> virtual get_matrix_elements(
        std::shared_ptr<const Derived> other, OperatorType type, int q = 0) const = 0;

    class Iterator {
    public:
        Iterator(typename ketvec_t::const_iterator it);
        bool operator!=(const Iterator &other) const;
        std::shared_ptr<const ket_t> operator*() const;
        Iterator &operator++();

    private:
        typename ketvec_t::const_iterator it;
    };

    Iterator begin() const;
    Iterator end() const;

    const Transformation<scalar_t> &get_transformation() const override;
    Transformation<scalar_t> get_rotator(real_t alpha, real_t beta, real_t gamma) const override;
    Sorting get_sorter(const std::vector<TransformationType> &labels) const override;
    std::vector<IndicesOfBlock>
    get_indices_of_blocks(const std::vector<TransformationType> &labels) const override;

    void perform_sorter_checks(const std::vector<TransformationType> &labels) const;
    void perform_blocks_checks(const std::set<TransformationType> &unique_labels) const;
    void get_sorter_without_checks(const std::vector<TransformationType> &labels,
                                   Sorting &transformation) const;
    void get_indices_of_blocks_without_checks(const std::set<TransformationType> &unique_labels,
                                              IndicesOfBlocksCreator &blocks) const;

    std::shared_ptr<const Derived>
    transformed(const Transformation<scalar_t> &transformation) const;
    std::shared_ptr<const Derived> transformed(const Sorting &transformation) const;

protected:
    Basis(ketvec_t &&kets);
    int get_ket_index_from_ket(std::shared_ptr<const ket_t> ket) const;
    ketvec_t kets;

private:
    const Derived &derived() const;

    struct hash {
        std::size_t operator()(const std::shared_ptr<const ket_t> &k) const;
    };

    struct equal_to {
        bool operator()(const std::shared_ptr<const ket_t> &lhs,
                        const std::shared_ptr<const ket_t> &rhs) const;
    };

    Transformation<scalar_t> coefficients;

    std::unordered_map<std::shared_ptr<const ket_t>, size_t, hash, equal_to> ket_to_ket_index;
    std::vector<size_t> ket_index_to_state_index;

    std::vector<real_t> state_index_to_quantum_number_f;
    std::vector<real_t> state_index_to_quantum_number_m;
    std::vector<Parity> state_index_to_parity;
    std::vector<size_t> state_index_to_ket_index;

    bool _has_quantum_number_f{true};
    bool _has_quantum_number_m{true};
    bool _has_parity{true};
};
} // namespace pairinteraction
