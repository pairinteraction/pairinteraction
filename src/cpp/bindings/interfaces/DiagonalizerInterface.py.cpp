#include "./DiagonalizerInterface.py.hpp"

#include "pairinteraction/interfaces/DiagonalizerInterface.hpp"

#include <nanobind/eigen/sparse.h>
#include <nanobind/nanobind.h>

namespace nb = nanobind;
using namespace pairinteraction;

template <typename T>
static void declare_diagonalizer_interface(nb::module_ &m, std::string const &type_name) {
    std::string pylass_name = "DiagonalizerInterface" + type_name;
    using real_t = typename DiagonalizerInterface<T>::real_t;
    nb::class_<DiagonalizerInterface<T>> pyclass(m, pylass_name.c_str());
    pyclass
        .def("eigh",
             nb::overload_cast<const Eigen::SparseMatrix<T, Eigen::RowMajor> &, int>(
                 &DiagonalizerInterface<T>::eigh, nb::const_))
        .def(
            "eigh",
            nb::overload_cast<const Eigen::SparseMatrix<T, Eigen::RowMajor> &, real_t, real_t, int>(
                &DiagonalizerInterface<T>::eigh, nb::const_));
}

void bind_diagonalizer_interface(nb::module_ &m) {
    declare_diagonalizer_interface<float>(m, "Float");
    declare_diagonalizer_interface<double>(m, "Double");
    declare_diagonalizer_interface<std::complex<float>>(m, "ComplexFloat");
    declare_diagonalizer_interface<std::complex<double>>(m, "ComplexDouble");
}
