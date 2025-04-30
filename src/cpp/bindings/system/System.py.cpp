// SPDX-FileCopyrightText: 2024 Pairinteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "./System.py.hpp"

#include "pairinteraction/basis/BasisAtom.hpp"
#include "pairinteraction/basis/BasisPair.hpp"
#include "pairinteraction/interfaces/DiagonalizerInterface.hpp"
#include "pairinteraction/system/GreenTensor.hpp"
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
#include <nanobind/stl/variant.h>
#include <nanobind/stl/vector.h>

namespace nb = nanobind;
using namespace nb::literals;
using namespace pairinteraction;

template <typename T>
static void declare_system(nb::module_ &m, std::string const &type_name) {
    using S = System<T>;
    using scalar_t = typename System<T>::scalar_t;

    std::string pyclass_name = "System" + type_name;

    nb::class_<System<T>, TransformationBuilderInterface<scalar_t>> pyclass(m,
                                                                            pyclass_name.c_str());
    pyclass.def("get_basis", &S::get_basis)
        .def("get_eigenbasis", &S::get_eigenbasis)
        .def("get_eigenenergies", &S::get_eigenenergies)
        .def("get_matrix", &S::get_matrix)
        .def("get_transformation", &S::get_transformation)
        .def("get_rotator", &S::get_rotator)
        .def("get_sorter", &S::get_sorter)
        .def("get_indices_of_blocks", &S::get_indices_of_blocks)
        .def("transform", nb::overload_cast<const Transformation<scalar_t> &>(&S::transform))
        .def("transform", nb::overload_cast<const Sorting &>(&S::transform))
        .def("diagonalize", &S::diagonalize, "diagonalizer"_a, "min_eigenenergy"_a = nb::none(),
             "max_eigenenergy"_a = nb::none(), "rtol"_a = 1e-6)
        .def("is_diagonal", &S::is_diagonal);
}

template <typename T>
static void declare_system_atom(nb::module_ &m, std::string const &type_name) {
    using S = SystemAtom<T>;
    using basis_t = typename SystemAtom<T>::basis_t;

    std::string pyclass_name = "SystemAtom" + type_name;

    nb::class_<S, System<S>> pyclass(m, pyclass_name.c_str());
    pyclass.def(nb::init<std::shared_ptr<const basis_t>>())
        .def("set_electric_field", &S::set_electric_field)
        .def("set_magnetic_field", &S::set_magnetic_field)
        .def("set_diamagnetism_enabled", &S::set_diamagnetism_enabled)
        .def("set_ion_distance_vector", &S::set_ion_distance_vector)
        .def("set_ion_charge", &S::set_ion_charge)
        .def("set_ion_interaction_order", &S::set_ion_interaction_order);
}

template <typename T>
static void declare_system_pair(nb::module_ &m, std::string const &type_name) {
    using S = SystemPair<T>;
    using basis_t = typename SystemPair<T>::basis_t;

    std::string pyclass_name = "SystemPair" + type_name;

    nb::class_<S, System<S>> pyclass(m, pyclass_name.c_str());
    pyclass.def(nb::init<std::shared_ptr<const basis_t>>())
        .def("set_interaction_order", &S::set_interaction_order)
        .def("set_distance_vector", &S::set_distance_vector)
        .def("set_green_tensor", &S::set_green_tensor);
}

template <typename T>
static void declare_green_tensor(nb::module_ &m, std::string const &type_name) {
    using CE = typename GreenTensor<T>::ConstantEntry;
    using OE = typename GreenTensor<T>::OmegaDependentEntry;
    using GT = GreenTensor<T>;

    std::string ce_name = "ConstantEntry" + type_name;
    std::string oe_name = "OmegaDependentEntry" + type_name;
    std::string gt_name = "GreenTensor" + type_name;

    nb::class_<CE>(m, ce_name.c_str())
        .def("row", &CE::row)
        .def("col", &CE::col)
        .def("val", &CE::val);

    nb::class_<OE>(m, oe_name.c_str())
        .def("row", &OE::row)
        .def("col", &OE::col)
        .def("val", &OE::val);

    nb::class_<GT>(m, gt_name.c_str())
        .def(nb::init<>())
        .def("set_entries",
             nb::overload_cast<int, int, const Eigen::MatrixX<T> &>(&GT::set_entries))
        .def("set_entries",
             nb::overload_cast<int, int, const std::vector<Eigen::MatrixX<T>> &,
                               const std::vector<double> &>(&GT::set_entries))
        .def("get_entries", nb::overload_cast<int, int>(&GT::get_entries, nb::const_));
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

    declare_green_tensor<double>(m, "Real");
    declare_green_tensor<std::complex<double>>(m, "Complex");
}
