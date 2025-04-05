// SPDX-FileCopyrightText: 2024 Pairinteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "./System.py.hpp"

#include "pairinteraction/basis/BasisAtom.hpp"
#include "pairinteraction/basis/BasisPair.hpp"
#include "pairinteraction/interfaces/DiagonalizerInterface.hpp"
#include "pairinteraction/system/System.hpp"
#include "pairinteraction/system/SystemAtom.hpp"
#include "pairinteraction/system/SystemPair.hpp"

#include <nanobind/eigen/dense.h>
#include <nanobind/eigen/sparse.h>
#include <nanobind/nanobind.h>
#include <nanobind/stl/array.h>
#include <nanobind/stl/complex.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/shared_ptr.h>
#include <nanobind/stl/vector.h>

namespace nb = nanobind;
using namespace nb::literals;
using namespace pairinteraction;

template <typename T>
static void declare_system(nb::module_ &m, std::string const &type_name) {
    std::string pyclass_name = "System" + type_name;
    using scalar_t = typename System<T>::scalar_t;
    nb::class_<System<T>, TransformationBuilderInterface<scalar_t>> pyclass(m,
                                                                            pyclass_name.c_str());
    pyclass.def("get_basis", &System<T>::get_basis)
        .def("get_eigenbasis", &System<T>::get_eigenbasis)
        .def("get_eigenenergies", &System<T>::get_eigenenergies)
        .def("get_matrix", &System<T>::get_matrix)
        .def("get_transformation", &System<T>::get_transformation)
        .def("get_rotator", &System<T>::get_rotator)
        .def("get_sorter", &System<T>::get_sorter)
        .def("get_indices_of_blocks", &System<T>::get_indices_of_blocks)
        .def("transform",
             nb::overload_cast<const Transformation<scalar_t> &>(&System<T>::transform))
        .def("transform", nb::overload_cast<const Sorting &>(&System<T>::transform))
        .def("diagonalize", &System<T>::diagonalize, "diagonalizer"_a,
             "min_eigenenergy"_a = nb::none(), "max_eigenenergy"_a = nb::none(), "rtol"_a = 1e-6)
        .def("is_diagonal", &System<T>::is_diagonal);
}

template <typename T>
static void declare_system_atom(nb::module_ &m, std::string const &type_name) {
    std::string pyclass_name = "SystemAtom" + type_name;
    using basis_t = typename SystemAtom<T>::basis_t;
    nb::class_<SystemAtom<T>, System<SystemAtom<T>>> pyclass(m, pyclass_name.c_str());
    pyclass.def(nb::init<std::shared_ptr<const basis_t>>())
        .def("set_electric_field", &SystemAtom<T>::set_electric_field)
        .def("set_magnetic_field", &SystemAtom<T>::set_magnetic_field)
        .def("set_diamagnetism_enabled", &SystemAtom<T>::set_diamagnetism_enabled)
        .def("set_ion_distance_vector", &SystemAtom<T>::set_ion_distance_vector)
        .def("set_ion_charge", &SystemAtom<T>::set_ion_charge)
        .def("set_ion_interaction_order", &SystemAtom<T>::set_ion_interaction_order);
}

template <typename T>
static void declare_system_pair(nb::module_ &m, std::string const &type_name) {
    std::string pyclass_name = "SystemPair" + type_name;
    using basis_t = typename SystemPair<T>::basis_t;
    nb::class_<SystemPair<T>, System<SystemPair<T>>> pyclass(m, pyclass_name.c_str());
    pyclass.def(nb::init<std::shared_ptr<const basis_t>>())
        .def("set_interaction_order", &SystemPair<T>::set_interaction_order)
        .def("set_distance_vector", &SystemPair<T>::set_distance_vector);
}

void bind_system(nb::module_ &m) {
    declare_system<SystemAtom<double>>(m, "SystemAtomReal");
    declare_system<SystemAtom<std::complex<double>>>(m, "SystemAtomComplex");
    declare_system_atom<double>(m, "Real");
    declare_system_atom<std::complex<double>>(m, "Complex");

    declare_system<SystemPair<double>>(m, "SystemPairReal");
    declare_system<SystemPair<std::complex<double>>>(m, "SystemPairComplex");
    declare_system_pair<double>(m, "Real");
    declare_system_pair<std::complex<double>>(m, "Complex");
}
