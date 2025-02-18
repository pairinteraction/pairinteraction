#include "./Diagonalizer.py.hpp"

#include "pairinteraction/diagonalizer/DiagonalizerEigen.hpp"
#include "pairinteraction/diagonalizer/DiagonalizerFeast.hpp"
#include "pairinteraction/diagonalizer/DiagonalizerLapacke.hpp"
#include "pairinteraction/diagonalizer/diagonalize.hpp"
#include "pairinteraction/enums/FloatType.hpp"
#include "pairinteraction/system/SystemAtom.hpp"
#include "pairinteraction/system/SystemPair.hpp"

#include <nanobind/eigen/sparse.h>
#include <nanobind/nanobind.h>
#include <nanobind/stl/complex.h>
#include <nanobind/stl/optional.h>

namespace nb = nanobind;
using namespace nb::literals;
using namespace pairinteraction;

template <typename T>
static void declare_diagonalizer_eigen(nb::module_ &m, std::string const &type_name) {
    std::string pyclass_name = "DiagonalizerEigen" + type_name;
    nb::class_<DiagonalizerEigen<T>, DiagonalizerInterface<T>> pyclass(m, pyclass_name.c_str());
    pyclass.def(nb::init<FloatType>(), "float_type"_a = FloatType::FLOAT64)
        .def("eigh",
             nb::overload_cast<const Eigen::SparseMatrix<T, Eigen::RowMajor> &, double>(
                 &DiagonalizerEigen<T>::eigh, nb::const_));
}

template <typename T>
static void declare_diagonalizer_feast(nb::module_ &m, std::string const &type_name) {
    std::string pyclass_name = "DiagonalizerFeast" + type_name;
    using real_t = typename DiagonalizerFeast<T>::real_t;
    nb::class_<DiagonalizerFeast<T>, DiagonalizerInterface<T>> pyclass(m, pyclass_name.c_str());
    pyclass.def(nb::init<int, FloatType>(), "m0"_a, "float_type"_a = FloatType::FLOAT64)
        .def("eigh",
             nb::overload_cast<const Eigen::SparseMatrix<T, Eigen::RowMajor> &, double>(
                 &DiagonalizerFeast<T>::eigh, nb::const_))
        .def("eigh",
             nb::overload_cast<const Eigen::SparseMatrix<T, Eigen::RowMajor> &,
                               std::optional<real_t>, std::optional<real_t>, double>(
                 &DiagonalizerFeast<T>::eigh, nb::const_));
}

template <typename T>
static void declare_diagonalizer_lapacke(nb::module_ &m, std::string const &type_name) {
    std::string pyclass_name = "DiagonalizerLapacke" + type_name;
    nb::class_<DiagonalizerLapacke<T>, DiagonalizerInterface<T>> pyclass(m, pyclass_name.c_str());
    pyclass.def(nb::init<FloatType>(), "float_type"_a = FloatType::FLOAT64)
        .def("eigh",
             nb::overload_cast<const Eigen::SparseMatrix<T, Eigen::RowMajor> &, double>(
                 &DiagonalizerLapacke<T>::eigh, nb::const_));
}

template <typename T>
static void declare_diagonalize(nb::module_ &m, std::string const &type_name) {
    std::string pyclass_name = "diagonalize" + type_name;
    using real_t = typename T::real_t;
    using scalar_t = typename T::scalar_t;
    m.def(
        pyclass_name.c_str(),
        [](nb::list pylist, // NOLINT
           const DiagonalizerInterface<scalar_t> &diagonalizer,
           std::optional<real_t> min_eigenvalue, std::optional<real_t> max_eigenvalue,
           double atol) {
            std::vector<T> systems;
            systems.reserve(pylist.size());
            for (auto h : pylist) {
                systems.push_back(nb::cast<T>(h));
            }
            diagonalize(systems, diagonalizer, min_eigenvalue, max_eigenvalue, atol);
            for (size_t i = 0; i < systems.size(); ++i) {
                pylist[i] = nb::cast(systems[i]);
            }
        },
        "systems"_a, "diagonalizer"_a, "min_eigenvalue"_a = nb::none(),
        "max_eigenvalue"_a = nb::none(), "atol"_a = 1e-6);
}

void bind_diagonalizer(nb::module_ &m) {
    declare_diagonalizer_eigen<double>(m, "Real");
    declare_diagonalizer_eigen<std::complex<double>>(m, "Complex");
    declare_diagonalizer_feast<double>(m, "Real");
    declare_diagonalizer_feast<std::complex<double>>(m, "Complex");
    declare_diagonalizer_lapacke<double>(m, "Real");
    declare_diagonalizer_lapacke<std::complex<double>>(m, "Complex");

    declare_diagonalize<SystemAtom<double>>(m, "SystemAtomReal");
    declare_diagonalize<SystemAtom<std::complex<double>>>(m, "SystemAtomComplex");

    declare_diagonalize<SystemPair<double>>(m, "SystemPairReal");
    declare_diagonalize<SystemPair<std::complex<double>>>(m, "SystemPairComplex");
}
