// SPDX-FileCopyrightText: 2024 PairInteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include "pairinteraction/interfaces/SorterBuilderInterface.hpp"
#include "pairinteraction/utils/eigen_assertion.hpp"
#include "pairinteraction/utils/eigen_compat.hpp"
#include "pairinteraction/utils/traits.hpp"

#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <memory>
#include <set>
#include <vector>

namespace pairinteraction {
enum class Parity : int;
enum class SorterType : unsigned char;

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
class Basis : public SorterBuilderInterface {
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
    const Eigen::SparseMatrix<scalar_t, Eigen::RowMajor> &get_coefficients() const;
    std::shared_ptr<const Derived>
    copy_with_coefficients(const Eigen::SparseMatrix<scalar_t, Eigen::RowMajor> &values) const;

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

    Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic>
    get_sorter(const std::vector<SorterType> &labels) const override;
    std::vector<IndicesOfBlock>
    get_indices_of_blocks(const std::vector<SorterType> &labels) const override;

    void perform_blocks_checks(const std::set<SorterType> &unique_labels) const;
    void get_sorter_without_checks(
        const std::vector<SorterType> &labels,
        Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> &sorter) const;
    void get_indices_of_blocks_without_checks(const std::set<SorterType> &unique_labels,
                                              IndicesOfBlocksCreator &blocks) const;

    std::shared_ptr<const Derived> canonicalized() const;
    bool is_canonical() const;
    virtual std::shared_ptr<const Derived> merge(std::shared_ptr<const Derived> other) const = 0;
    std::shared_ptr<const Derived>
    transformed(const Eigen::SparseMatrix<scalar_t, Eigen::RowMajor> &transformation) const;
    std::shared_ptr<const Derived>
    transformed(const Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> &sorter) const;

protected:
    Basis(ketvec_t &&kets);
    ketvec_t kets;

private:
    const Derived &derived() const;

    Eigen::SparseMatrix<scalar_t, Eigen::RowMajor> coefficients;

    std::vector<real_t> state_index_to_quantum_number_f;
    std::vector<real_t> state_index_to_quantum_number_m;
    std::vector<Parity> state_index_to_parity;
};
} // namespace pairinteraction
