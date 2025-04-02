// SPDX-FileCopyrightText: 2024 Pairinteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "./TransformationBuilderInterface.py.hpp"

#include "pairinteraction/interfaces/TransformationBuilderInterface.hpp"

#include <nanobind/eigen/dense.h>
#include <nanobind/eigen/sparse.h>
#include <nanobind/nanobind.h>
#include <nanobind/stl/array.h>
#include <nanobind/stl/complex.h>
#include <nanobind/stl/vector.h>

namespace nb = nanobind;
using namespace nb::literals;
using namespace pairinteraction;

template <typename T>
static void declare_transformation(nb::module_ &m, std::string const &type_name) {
    std::string pyclass_name = "Transformation" + type_name;
    nb::class_<Transformation<T>> pyclass(m, pyclass_name.c_str());
    pyclass.def(nb::init<>())
        .def_rw("matrix", &Transformation<T>::matrix)
        .def_rw("transformation_type", &Transformation<T>::transformation_type);
}

static void declare_sorting(nb::module_ &m) {
    nb::class_<Sorting> pyclass(m, "Sorting");
    pyclass.def(nb::init<>())
        .def_rw("matrix", &Sorting::matrix)
        .def_rw("transformation_type", &Sorting::transformation_type);
}

static void declare_indices_of_blocks(nb::module_ &m) {
    nb::class_<IndicesOfBlock> pyclass(m, "IndicesOfBlock");
    pyclass.def(nb::init<size_t, size_t>(), "start"_a, "end"_a)
        .def_rw("start", &IndicesOfBlock::start)
        .def_rw("end", &IndicesOfBlock::end);
}

static void declare_indices_of_blocks_creator(nb::module_ &m) {
    nb::class_<IndicesOfBlocksCreator> pyclass(m, "IndicesOfBlocksCreator");
}

template <typename T>
static void declare_transformation_builder_interface(nb::module_ &m, std::string const &type_name) {
    std::string pyclass_name = "TransformationBuilderInterface" + type_name;
    using real_t = typename TransformationBuilderInterface<T>::real_t;
    nb::class_<TransformationBuilderInterface<T>> pyclass(m, pyclass_name.c_str());
    pyclass.def("get_rotator",
                nb::overload_cast<const std::array<real_t, 3> &, const std::array<real_t, 3> &>(
                    &TransformationBuilderInterface<T>::get_rotator, nb::const_));
}

void bind_transformation_builder_interface(nb::module_ &m) {
    declare_transformation<double>(m, "Real");
    declare_transformation<std::complex<double>>(m, "Complex");
    declare_sorting(m);
    declare_indices_of_blocks(m);
    declare_indices_of_blocks_creator(m);
    declare_transformation_builder_interface<double>(m, "Real");
    declare_transformation_builder_interface<std::complex<double>>(m, "Complex");
}
