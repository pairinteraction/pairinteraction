#include "./Diagonalizer.py.hpp"

#include "pairinteraction/diagonalizer/DiagonalizerEigen.hpp"
#include "pairinteraction/diagonalizer/DiagonalizerFeast.hpp"
#include "pairinteraction/diagonalizer/DiagonalizerLapacke.hpp"
#include "pairinteraction/diagonalizer/diagonalize.hpp"
#include "pairinteraction/system/SystemAtom.hpp"

#include <nanobind/eigen/sparse.h>
#include <nanobind/nanobind.h>

namespace nb = nanobind;
using namespace nb::literals;
using namespace pairinteraction;

template <typename T>
static void declare_diagonalizer_eigen(nb::module_ &m, std::string const &type_name) {
    std::string pyclass_name = "DiagonalizerEigen" + type_name;
    using real_t = typename DiagonalizerEigen<T>::real_t;
    nb::class_<DiagonalizerEigen<T>, DiagonalizerInterface<T>> pyclass(m, pyclass_name.c_str());
    pyclass.def(nb::init<>())
        .def("eigh",
             nb::overload_cast<const Eigen::SparseMatrix<T, Eigen::RowMajor> &, int>(
                 &DiagonalizerEigen<T>::eigh, nb::const_));
}

template <typename T>
static void declare_diagonalizer_feast(nb::module_ &m, std::string const &type_name) {
    std::string pyclass_name = "DiagonalizerFeast" + type_name;
    using real_t = typename DiagonalizerFeast<T>::real_t;
    nb::class_<DiagonalizerFeast<T>, DiagonalizerInterface<T>> pyclass(m, pyclass_name.c_str());
    pyclass.def(nb::init<int>())
        .def("eigh",
             nb::overload_cast<const Eigen::SparseMatrix<T, Eigen::RowMajor> &, int>(
                 &DiagonalizerFeast<T>::eigh, nb::const_))
        .def(
            "eigh",
            nb::overload_cast<const Eigen::SparseMatrix<T, Eigen::RowMajor> &, real_t, real_t, int>(
                &DiagonalizerFeast<T>::eigh, nb::const_));
}

template <typename T>
static void declare_diagonalizer_lapacke(nb::module_ &m, std::string const &type_name) {
    std::string pyclass_name = "DiagonalizerLapacke" + type_name;
    using real_t = typename DiagonalizerLapacke<T>::real_t;
    nb::class_<DiagonalizerLapacke<T>, DiagonalizerInterface<T>> pyclass(m, pyclass_name.c_str());
    pyclass.def(nb::init<>())
        .def("eigh",
             nb::overload_cast<const Eigen::SparseMatrix<T, Eigen::RowMajor> &, int>(
                 &DiagonalizerLapacke<T>::eigh, nb::const_));
}

template <typename T>
static void declare_diagonalize(nb::module_ &m, std::string const & /* type_name */) {
    using real_t = typename T::real_t;
    using scalar_t = typename T::scalar_t;
    m.def(
        "diagonalize",
        [](nb::list pylist, // NOLINT
           const DiagonalizerInterface<scalar_t> &diagonalizer, int precision) {
            std::vector<std::reference_wrapper<T>> systems;
            systems.reserve(pylist.size());
            for (auto h : pylist) {
                T &obj = nb::cast<T &>(h);
                systems.push_back(std::ref(obj));
            }
            return diagonalize(systems, diagonalizer, precision);
        },
        "systems"_a, "diagonalizer"_a, "precision"_a = 12);
    m.def(
        "diagonalize",
        [](nb::list pylist, // NOLINT
           const DiagonalizerInterface<scalar_t> &diagonalizer, real_t min_eigenvalue,
           real_t max_eigenvalue, int precision) {
            std::vector<std::reference_wrapper<T>> systems;
            systems.reserve(pylist.size());
            for (auto h : pylist) {
                T &obj = nb::cast<T &>(h);
                systems.push_back(std::ref(obj));
            }
            return diagonalize(systems, diagonalizer, min_eigenvalue, max_eigenvalue, precision);
        },
        "systems"_a, "diagonalizer"_a, "min_eigenvalue"_a, "max_eigenvalue"_a, "precision"_a = 12);
}

void bind_diagonalizer(nb::module_ &m) {
    declare_diagonalizer_eigen<float>(m, "Float");
    declare_diagonalizer_eigen<double>(m, "Double");
    declare_diagonalizer_eigen<std::complex<float>>(m, "ComplexFloat");
    declare_diagonalizer_eigen<std::complex<double>>(m, "ComplexDouble");
    declare_diagonalizer_feast<float>(m, "Float");
    declare_diagonalizer_feast<double>(m, "Double");
    declare_diagonalizer_feast<std::complex<float>>(m, "ComplexFloat");
    declare_diagonalizer_feast<std::complex<double>>(m, "ComplexDouble");
    declare_diagonalizer_lapacke<float>(m, "Float");
    declare_diagonalizer_lapacke<double>(m, "Double");
    declare_diagonalizer_lapacke<std::complex<float>>(m, "ComplexFloat");
    declare_diagonalizer_lapacke<std::complex<double>>(m, "ComplexDouble");
    declare_diagonalize<SystemAtom<float>>(m, "SystemAtomFloat");
    declare_diagonalize<SystemAtom<double>>(m, "SystemAtomDouble");
    declare_diagonalize<SystemAtom<std::complex<float>>>(m, "SystemAtomComplexFloat");
    declare_diagonalize<SystemAtom<std::complex<double>>>(m, "SystemAtomComplexDouble");
}
