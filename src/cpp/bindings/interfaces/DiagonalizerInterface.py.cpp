// SPDX-FileCopyrightText: 2024 Pairinteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "./DiagonalizerInterface.py.hpp"

#include "pairinteraction/interfaces/DiagonalizerInterface.hpp"

#include <nanobind/eigen/sparse.h>
#include <nanobind/nanobind.h>
#include <nanobind/stl/complex.h>
#include <nanobind/stl/optional.h>

namespace nb = nanobind;
using namespace pairinteraction;

template <typename T>
static void declare_diagonalizer_interface(nb::module_ &m, std::string const &type_name) {
    std::string pylass_name = "DiagonalizerInterface" + type_name;
    using real_t = typename DiagonalizerInterface<T>::real_t;
    nb::class_<DiagonalizerInterface<T>> pyclass(m, pylass_name.c_str());
    pyclass
        .def("eigh",
             nb::overload_cast<const Eigen::SparseMatrix<T, Eigen::RowMajor> &, double>(
                 &DiagonalizerInterface<T>::eigh, nb::const_))
        .def("eigh",
             nb::overload_cast<const Eigen::SparseMatrix<T, Eigen::RowMajor> &,
                               std::optional<real_t>, std::optional<real_t>, double>(
                 &DiagonalizerInterface<T>::eigh, nb::const_));
}

template <typename T>
static void declare_eigen_system_h(nb::module_ &m, std::string const &type_name) {
    std::string pylass_name = "EigenSystemH" + type_name;
    nb::class_<EigenSystemH<T>> pyclass(m, pylass_name.c_str());
    pyclass.def_rw("eigenvectors", &EigenSystemH<T>::eigenvectors)
        .def_rw("eigenvalues", &EigenSystemH<T>::eigenvalues);
}

void bind_diagonalizer_interface(nb::module_ &m) {
    declare_diagonalizer_interface<double>(m, "Real");
    declare_diagonalizer_interface<std::complex<double>>(m, "Complex");
    declare_eigen_system_h<double>(m, "Real");
    declare_eigen_system_h<std::complex<double>>(m, "Complex");
}
