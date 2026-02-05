// SPDX-FileCopyrightText: 2024 PairInteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "./Diagonalizer.py.hpp"

#include "pairinteraction/diagonalize/DiagonalizerEigen.hpp"
#include "pairinteraction/diagonalize/DiagonalizerFeast.hpp"
#include "pairinteraction/diagonalize/DiagonalizerLapackeEvd.hpp"
#include "pairinteraction/diagonalize/DiagonalizerLapackeEvr.hpp"
#include "pairinteraction/diagonalize/DiagonalizerSpectra.hpp"
#include "pairinteraction/diagonalize/diagonalize.hpp"
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
        .def("eigh_full", &DiagonalizerEigen<T>::eigh_full);
}

template <typename T>
static void declare_diagonalizer_lapacke_evd(nb::module_ &m, std::string const &type_name) {
    std::string pyclass_name = "DiagonalizerLapackeEvd" + type_name;
    nb::class_<DiagonalizerLapackeEvd<T>, DiagonalizerInterface<T>> pyclass(m,
                                                                            pyclass_name.c_str());
    pyclass.def(nb::init<FloatType>(), "float_type"_a = FloatType::FLOAT64)
        .def("eigh_full", &DiagonalizerLapackeEvd<T>::eigh_full);
}

template <typename T>
static void declare_diagonalizer_lapacke_evr(nb::module_ &m, std::string const &type_name) {
    std::string pyclass_name = "DiagonalizerLapackeEvr" + type_name;
    nb::class_<DiagonalizerLapackeEvr<T>, DiagonalizerInterface<T>> pyclass(m,
                                                                            pyclass_name.c_str());
    pyclass.def(nb::init<FloatType>(), "float_type"_a = FloatType::FLOAT64)
        .def("eigh_full", &DiagonalizerLapackeEvr<T>::eigh_full)
        .def("eigh_range", &DiagonalizerLapackeEvr<T>::eigh_range);
}

template <typename T>
static void declare_diagonalizer_feast(nb::module_ &m, std::string const &type_name) {
    std::string pyclass_name = "DiagonalizerFeast" + type_name;
    using real_t = typename DiagonalizerFeast<T>::real_t;
    nb::class_<DiagonalizerFeast<T>, DiagonalizerInterface<T>> pyclass(m, pyclass_name.c_str());
    pyclass.def(nb::init<int, FloatType>(), "m0"_a, "float_type"_a = FloatType::FLOAT64)
        .def("eigh_full", &DiagonalizerFeast<T>::eigh_full)
        .def("eigh_range", &DiagonalizerFeast<T>::eigh_range);
}

template <typename T>
static void declare_diagonalizer_spectra(nb::module_ &m, std::string const &type_name) {
    std::string pyclass_name = "DiagonalizerSpectra" + type_name;
    using real_t = typename DiagonalizerSpectra<T>::real_t;
    nb::class_<DiagonalizerSpectra<T>, DiagonalizerInterface<T>> pyclass(m, pyclass_name.c_str());
    pyclass
        .def(nb::init<std::optional<real_t>, FloatType>(), "ncv"_a = nb::none(),
             "float_type"_a = FloatType::FLOAT64)
        .def("eigh_full", &DiagonalizerSpectra<T>::eigh_full)
        .def("eigh_range", &DiagonalizerSpectra<T>::eigh_range)
        .def("eigh_shift_invert", &DiagonalizerSpectra<T>::eigh_shift_invert);
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
           double rtol) {
            std::vector<std::reference_wrapper<T>> systems;
            systems.reserve(pylist.size());
            for (nb::handle_t<T> &&h : pylist) {
                systems.push_back(nb::cast<T &>(h));
            }
            diagonalize(systems, diagonalizer, min_eigenvalue, max_eigenvalue, rtol);
        },
        "systems"_a, "diagonalizer"_a, "min_eigenvalue"_a = nb::none(),
        "max_eigenvalue"_a = nb::none(), "rtol"_a = 1e-6);
}

void bind_diagonalizer(nb::module_ &m) {
    declare_diagonalizer_eigen<double>(m, "Real");
    declare_diagonalizer_eigen<std::complex<double>>(m, "Complex");
    declare_diagonalizer_feast<double>(m, "Real");
    declare_diagonalizer_feast<std::complex<double>>(m, "Complex");
    declare_diagonalizer_lapacke_evd<double>(m, "Real");
    declare_diagonalizer_lapacke_evd<std::complex<double>>(m, "Complex");
    declare_diagonalizer_lapacke_evr<double>(m, "Real");
    declare_diagonalizer_lapacke_evr<std::complex<double>>(m, "Complex");
    declare_diagonalizer_spectra<double>(m, "Real");
    declare_diagonalizer_spectra<std::complex<double>>(m, "Complex");

    declare_diagonalize<SystemAtom<double>>(m, "SystemAtomReal");
    declare_diagonalize<SystemAtom<std::complex<double>>>(m, "SystemAtomComplex");

    declare_diagonalize<SystemPair<double>>(m, "SystemPairReal");
    declare_diagonalize<SystemPair<std::complex<double>>>(m, "SystemPairComplex");
}
