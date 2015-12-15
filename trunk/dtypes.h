#ifndef DTYPES_H
#define DTYPES_H

#include "utils.hpp"

#include <vector>
#include <array>
#include <iostream>
#include <iomanip>
#include <Eigen/Sparse>
#include <Eigen/Dense>

#define REAL_T MPI_FLOAT
#define IDX_T MPI_UNSIGNED
#define BYTE_T MPI_UNSIGNED_CHAR
#define ORDER Eigen::ColMajor

typedef float real_t;
typedef std::complex<real_t> complex_t;
typedef uint32_t idx_t;
typedef float storage_real_t;
typedef int32_t storage_idx_t;
typedef int32_t eigen_idx_t;

typedef uint8_t byte_t;
typedef std::vector<byte_t> bytes_t;
typedef std::nullptr_t invalid_t;

#include "config.h"

typedef Eigen::Triplet<scalar_t> eigen_triplet_t;
typedef Eigen::Triplet<real_t> eigen_triplet_real_t;
typedef Eigen::SparseMatrix<scalar_t, ORDER, eigen_idx_t> eigen_sparse_t;
typedef Eigen::SparseMatrix<real_t, ORDER, eigen_idx_t> eigen_sparse_real_t;
typedef Eigen::SparseMatrix<scalar_t, ORDER, eigen_idx_t>::InnerIterator eigen_iterator_t;
typedef Eigen::SparseMatrix<real_t, ORDER, eigen_idx_t>::InnerIterator eigen_iterator_real_t;
typedef Eigen::Matrix<scalar_t,Eigen::Dynamic,Eigen::Dynamic,ORDER> eigen_dense_t;
typedef Eigen::Matrix<scalar_t,Eigen::Dynamic,1,ORDER> eigen_vector_t;
typedef Eigen::Matrix<real_t,Eigen::Dynamic,1,ORDER> eigen_vector_real_t;


class Triple {
public:
    Triple() : row(0), col(0), val(0) { }
    Triple(idx_t row, idx_t col, scalar_t val) : row(row), col(col), val(val) { }
    idx_t row;
    idx_t col;
    scalar_t val;
};

class State {
public:
    State(idx_t idx) : idx(idx) { }
    idx_t idx;
};

class StateOne : public State {
public:
    StateOne() : State(0), n(0), l(0), s(0), j(0), m(0) { }
    StateOne(idx_t idx, int n, int l, float s, float j, float m) : State(idx), n(n), l(l), s(s), j(j), m(m) { }
    StateOne(int n, int l, float s, float j, float m) : State(0), n(n), l(l), s(s), j(j), m(m) { }
    friend std::ostream& operator<< (std::ostream &out, const StateOne &state) {
        out << "i  =" << std::setw(5) << state.idx << ",   ";
        out << "n  =" << std::setw(3) << state.n << ",   ";
        out << "l  =" << std::setw(2) << state.l << ",   ";
        out << "s  =" << std::setprecision(2) << std::setw(4) << state.s << ",   ";
        out << "j  =" << std::setprecision(2) << std::setw(4) << state.j << ",   ";
        out << "m  =" << std::setprecision(2) << std::setw(4) << state.m;
        return out;
    }

    friend bool operator== (const StateOne& s1, const StateOne& s2)
    {
        return ((s1.n == s2.n) && (s1.l == s2.l)  && (s1.s == s2.s)  && (s1.j == s2.j)  && (s1.m == s2.m));
    }

    int n, l;
    float s, j, m;
};

class StateTwo : public State {
public:
    StateTwo() : State(0), n({{0,0}}), l({{0,0}}), s({{0,0}}), j({{0,0}}), m({{0,0}}) { }
    StateTwo(idx_t idx, std::array<int, 2> n, std::array<int, 2> l, std::array<float, 2> s, std::array<float, 2> j, std::array<float, 2> m) : State(idx), n(n), l(l), s(s), j(j), m(m) { }
    StateTwo(std::array<int, 2> n, std::array<int, 2> l, std::array<float, 2> s, std::array<float, 2> j, std::array<float, 2> m) : State(0), n(n), l(l), s(s), j(j), m(m) { }
    StateTwo(idx_t idx, const StateOne &a, const StateOne &b) : State(idx), n({{a.n,b.n}}), l({{a.l,b.l}}), s({{a.s,b.s}}), j({{a.j,b.j}}), m({{a.m,b.m}}) { }
    StateTwo(const StateOne &a, const StateOne &b) : State(0), n({{a.n,b.n}}), l({{a.l,b.l}}), s({{a.s,b.s}}), j({{a.j,b.j}}), m({{a.m,b.m}}) { }
    friend std::ostream& operator<< (std::ostream &out, const StateTwo &state) {
        out << "i  =" << std::setw(5) << state.idx << ",   ";
        for (size_t i = 0; i < 2; ++i) {
            out << "n" << i << " =" << std::setw(3) << state.n[i] << ",   ";
            out << "l" << i << " =" << std::setw(2) << state.l[i] << ",   ";
            out << "s" << i << " =" << std::setprecision(2) << std::setw(4) << state.s[i] << ",   ";
            out << "j" << i << " =" << std::setprecision(2) << std::setw(4) << state.j[i] << ",   ";
            out << "m" << i << " =" << std::setprecision(2) << std::setw(4) << state.m[i];
            if (i == 0) out << ",   ";
        }
        return out;
    }
    friend bool operator== (const StateTwo& s1, const StateTwo& s2)
    {
        return ((s1.n[0] == s2.n[0]) && (s1.l[0] == s2.l[0])  && (s1.s[0] == s2.s[0])  && (s1.j[0] == s2.j[0])  && (s1.m[0] == s2.m[0]) &&
                (s1.n[1] == s2.n[1]) && (s1.l[1] == s2.l[1])  && (s1.s[1] == s2.s[1])  && (s1.j[1] == s2.j[1])  && (s1.m[1] == s2.m[1]));
    }
    StateOne first() {
        return StateOne(idx, n[0], l[0], s[0], j[0], m[0]);
    }
    StateOne second() {
        return StateOne(idx, n[1], l[1], s[1], j[1], m[1]);
    }
    StateTwo order() {
        if ((n[0] < n[1]) || ((n[0] == n[1]) &&
           ((l[0] < l[1]) || ((l[0] == l[1]) &&
           ((s[0] < s[1]) || ((s[0] == s[1]) &&
           ((j[0] < j[1]) || ((j[0] == j[1]) &&
           (m[0] <= m[1]))))))))) {
            return *this;
        } else {
            return StateTwo(this->second(),this->first());
        }
    }
    std::array<int, 2> n, l;
    std::array<float, 2> s, j, m;
};

namespace std {
template <> struct hash<StateOne>
{
    size_t operator()(const StateOne & s) const
    {
        size_t seed = 0;
        utils::hash_combine(seed, s.n);
        utils::hash_combine(seed, s.l);
        utils::hash_combine(seed, s.s);
        utils::hash_combine(seed, s.j);
        utils::hash_combine(seed, s.m);
        return seed;
    }
};
template <> struct hash<StateTwo>
{
    size_t operator()(const StateTwo & s) const
    {
        size_t seed = 0;
        utils::hash_combine(seed, s.n[0]);
        utils::hash_combine(seed, s.l[0]);
        utils::hash_combine(seed, s.s[0]);
        utils::hash_combine(seed, s.j[0]);
        utils::hash_combine(seed, s.m[0]);
        utils::hash_combine(seed, s.n[1]);
        utils::hash_combine(seed, s.l[1]);
        utils::hash_combine(seed, s.s[1]);
        utils::hash_combine(seed, s.j[1]);
        utils::hash_combine(seed, s.m[1]);
        return seed;
    }
};
}

#endif
