/*
 * Copyright (c) 2016 Sebastian Weber, Henri Menke. All rights reserved.
 *
 * This file is part of the pairinteraction library.
 *
 * The pairinteraction library is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * The pairinteraction library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with the pairinteraction library. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef HAMILTONIANMATRIX_H
#define HAMILTONIANMATRIX_H

#include "Basisnames.hpp"
#include "Serializable.hpp"
#include "StateOld.hpp"
#include "dtypes.hpp"
#include "utils.hpp"

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <exception>
#include <memory>
#include <stdexcept>
#include <vector>

const uint8_t csr_not_csc = 0x01;      // xxx0: csc, xxx1: csr
const uint8_t complex_not_real = 0x02; // xx0x: real, xx1x: complex

template <typename Scalar>
class Hamiltonianmatrix;

template <typename Scalar>
Hamiltonianmatrix<Scalar> operator+(Hamiltonianmatrix<Scalar> lhs,
                                    const Hamiltonianmatrix<Scalar> &rhs);
template <typename Scalar>
Hamiltonianmatrix<Scalar> operator-(Hamiltonianmatrix<Scalar> lhs,
                                    const Hamiltonianmatrix<Scalar> &rhs);
template <typename Scalar, typename T>
Hamiltonianmatrix<Scalar> operator*(const T &lhs, Hamiltonianmatrix<Scalar> rhs);
template <typename Scalar, typename T>
Hamiltonianmatrix<Scalar> operator*(Hamiltonianmatrix<Scalar> lhs, const T &rhs);

template <typename Scalar>
Hamiltonianmatrix<Scalar>
combine(const Hamiltonianmatrix<Scalar> &lhs, const Hamiltonianmatrix<Scalar> &rhs,
        const double &deltaE, const std::shared_ptr<BasisnamesTwo> &basis_two, const Symmetry &sym);
template <typename Scalar>
void energycutoff(const Hamiltonianmatrix<Scalar> &lhs, const Hamiltonianmatrix<Scalar> &rhs,
                  const double &deltaE, std::vector<bool> &necessary);

template <typename Scalar>
class Hamiltonianmatrix : public Serializable {
public:
    Hamiltonianmatrix();
    Hamiltonianmatrix(const Eigen::SparseMatrix<Scalar> &entries,
                      const Eigen::SparseMatrix<Scalar> &basis);
    Hamiltonianmatrix(size_t szBasis, size_t szEntries);
    Eigen::SparseMatrix<Scalar> &entries();
    const Eigen::SparseMatrix<Scalar> &entries() const;
    Eigen::SparseMatrix<Scalar> &basis();
    const Eigen::SparseMatrix<Scalar> &basis() const;
    size_t num_basisvectors() const;
    size_t num_coordinates() const;
    void addBasis(idx_t row, idx_t col, Scalar val);
    void addEntries(idx_t row, idx_t col, Scalar val);
    void compress(size_t nBasis, size_t nCoordinates);
    std::vector<Hamiltonianmatrix> findSubs() const;
    Hamiltonianmatrix abs() const;
    Hamiltonianmatrix changeBasis(const Eigen::SparseMatrix<Scalar> &basis) const;
    void applyCutoff(double cutoff);
    void findUnnecessaryStates(std::vector<bool> &isNecessaryCoordinate) const;
    void removeUnnecessaryBasisvectors(const std::vector<bool> &isNecessaryCoordinate);
    void removeUnnecessaryBasisvectors();
    void removeUnnecessaryStates(const std::vector<bool> &isNecessaryCoordinate);
    Hamiltonianmatrix getBlock(const std::vector<ptrdiff_t> &indices);
    void diagonalize();
    friend Hamiltonianmatrix operator+<>(Hamiltonianmatrix lhs, const Hamiltonianmatrix &rhs);
    friend Hamiltonianmatrix operator-<>(Hamiltonianmatrix lhs, const Hamiltonianmatrix &rhs);
    template <typename S, typename T>
    friend Hamiltonianmatrix<S> operator*(const T &lhs, Hamiltonianmatrix<S> rhs); // NOLINT
    template <typename S, typename T>
    friend Hamiltonianmatrix<S> operator*(Hamiltonianmatrix<S> lhs, const T &rhs); // NOLINT
    Hamiltonianmatrix &operator+=(const Hamiltonianmatrix &rhs);
    Hamiltonianmatrix &operator-=(const Hamiltonianmatrix &rhs);
    bytes_t &serialize() override;
    void doSerialization();
    void deserialize(bytes_t &bytesin) override;
    void doDeserialization();
    uint64_t hashEntries();
    uint64_t hashBasis();
    void save(const std::string &fname);
    bool load(const std::string &fname);
    friend Hamiltonianmatrix combine<>(const Hamiltonianmatrix &lhs, const Hamiltonianmatrix &rhs,
                                       const double &deltaE,
                                       const std::shared_ptr<BasisnamesTwo> &basis_two,
                                       const Symmetry &sym);
    friend void energycutoff<>(const Hamiltonianmatrix &lhs, const Hamiltonianmatrix &rhs,
                               const double &deltaE, std::vector<bool> &necessary);

    template <typename T, typename std::enable_if<utils::is_complex<T>::value>::type * = nullptr>
    void mergeComplex(std::vector<storage_double> &real, std::vector<storage_double> &imag,
                      std::vector<T> &complex) {
        std::vector<storage_double>::iterator real_it, imag_it;
        complex.reserve(real.size());
        for (real_it = real.begin(), imag_it = imag.begin(); real_it != real.end();
             ++real_it, ++imag_it) {
            complex.push_back(T(*real_it, *imag_it));
        }
    }

    template <typename T, typename std::enable_if<!utils::is_complex<T>::value>::type * = nullptr>
    void mergeComplex(std::vector<storage_double> &real, std::vector<storage_double> &imag,
                      std::vector<T> &complex) {
        (void)imag;
        complex = real;
    }

    template <typename T, typename std::enable_if<utils::is_complex<T>::value>::type * = nullptr>
    void splitComplex(std::vector<storage_double> &real, std::vector<storage_double> &imag,
                      std::vector<T> &complex) {
        real.reserve(complex.size());
        imag.reserve(imag.size());
        for (auto complex_it = complex.begin(); complex_it != complex.end(); ++complex_it) {
            real.push_back(complex_it->real());
            imag.push_back(complex_it->imag());
        }
    }

    template <typename T, typename std::enable_if<!utils::is_complex<T>::value>::type * = nullptr>
    void splitComplex(std::vector<storage_double> &real, std::vector<storage_double> &imag,
                      std::vector<T> &complex) {
        imag = std::vector<storage_double>();
        real = complex; // std::vector<storage_double>(complex.begin(),complex.end());
    }

protected:
    Eigen::SparseMatrix<Scalar> entries_;
    Eigen::SparseMatrix<Scalar> basis_;
    bytes_t bytes;
    std::vector<Eigen::Triplet<Scalar>> triplets_basis;
    std::vector<Eigen::Triplet<Scalar>> triplets_entries;
};

extern template class Hamiltonianmatrix<std::complex<double>>;
extern template Hamiltonianmatrix<std::complex<double>>
operator+(Hamiltonianmatrix<std::complex<double>> lhs,
          const Hamiltonianmatrix<std::complex<double>> &rhs);
extern template Hamiltonianmatrix<std::complex<double>>
operator-(Hamiltonianmatrix<std::complex<double>> lhs,
          const Hamiltonianmatrix<std::complex<double>> &rhs);
extern template Hamiltonianmatrix<std::complex<double>>
operator*(const std::complex<double> &lhs, Hamiltonianmatrix<std::complex<double>> rhs);
extern template Hamiltonianmatrix<std::complex<double>>
operator*(Hamiltonianmatrix<std::complex<double>> lhs, const std::complex<double> &rhs);
extern template Hamiltonianmatrix<std::complex<double>>
operator*(const double &lhs, Hamiltonianmatrix<std::complex<double>> rhs);
extern template Hamiltonianmatrix<std::complex<double>>
operator*(Hamiltonianmatrix<std::complex<double>> lhs, const double &rhs);
extern template Hamiltonianmatrix<std::complex<double>>
combine(const Hamiltonianmatrix<std::complex<double>> &lhs,
        const Hamiltonianmatrix<std::complex<double>> &rhs, const double &deltaE,
        const std::shared_ptr<BasisnamesTwo> &basis_two, const Symmetry &sym);
extern template void energycutoff(const Hamiltonianmatrix<std::complex<double>> &lhs,
                                  const Hamiltonianmatrix<std::complex<double>> &rhs,
                                  const double &deltaE, std::vector<bool> &necessary);
extern template class Hamiltonianmatrix<double>;
extern template Hamiltonianmatrix<double> operator+(Hamiltonianmatrix<double> lhs,
                                                    const Hamiltonianmatrix<double> &rhs);
extern template Hamiltonianmatrix<double> operator-(Hamiltonianmatrix<double> lhs,
                                                    const Hamiltonianmatrix<double> &rhs);
extern template Hamiltonianmatrix<double> operator*(const double &lhs,
                                                    Hamiltonianmatrix<double> rhs);
extern template Hamiltonianmatrix<double> operator*(Hamiltonianmatrix<double> lhs,
                                                    const double &rhs);
extern template Hamiltonianmatrix<double>
combine(const Hamiltonianmatrix<double> &lhs, const Hamiltonianmatrix<double> &rhs,
        const double &deltaE, const std::shared_ptr<BasisnamesTwo> &basis_two, const Symmetry &sym);
extern template void energycutoff(const Hamiltonianmatrix<double> &lhs,
                                  const Hamiltonianmatrix<double> &rhs, const double &deltaE,
                                  std::vector<bool> &necessary);

#endif // HAMILTONIANMATRIX_H
