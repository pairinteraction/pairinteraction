// SPDX-FileCopyrightText: 2024 Pairinteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "./Operator.py.hpp"

#include "pairinteraction/basis/BasisAtom.hpp"
#include "pairinteraction/basis/BasisPair.hpp"
#include "pairinteraction/operator/Operator.hpp"
#include "pairinteraction/operator/OperatorAtom.hpp"
#include "pairinteraction/operator/OperatorPair.hpp"

#include <nanobind/eigen/sparse.h>
#include <nanobind/nanobind.h>
#include <nanobind/operators.h>
#include <nanobind/stl/complex.h>
#include <nanobind/stl/shared_ptr.h>
#include <nanobind/stl/vector.h>

namespace nb = nanobind;
using namespace pairinteraction;

template <typename T>
static void declare_operator(nb::module_ &m, std::string const &type_name) {
    std::string pyclass_name = "Operator" + type_name;
    using basis_t = typename Operator<T>::basis_t;
    using scalar_t = typename Operator<T>::scalar_t;
    nb::class_<Operator<T>, TransformationBuilderInterface<scalar_t>> pyclass(m,
                                                                              pyclass_name.c_str());
    pyclass.def(nb::init<std::shared_ptr<const basis_t>>())
        .def("get_basis", nb::overload_cast<>(&Operator<T>::get_basis, nb::const_))
        .def("get_matrix", nb::overload_cast<>(&Operator<T>::get_matrix, nb::const_))
        .def("get_transformation", &Operator<T>::get_transformation)
        .def("get_rotator", &Operator<T>::get_rotator)
        .def("get_sorter", &Operator<T>::get_sorter)
        .def("get_indices_of_blocks", &Operator<T>::get_indices_of_blocks)
        .def("transformed",
             nb::overload_cast<const Transformation<scalar_t> &>(&Operator<T>::transformed,
                                                                 nb::const_))
        .def("transformed",
             nb::overload_cast<const Sorting &>(&Operator<T>::transformed, nb::const_))
        .def(scalar_t() * nb::self)
        .def(nb::self * scalar_t())
        .def(nb::self / scalar_t())
        .def(nb::self + nb::self)
        .def(nb::self - nb::self); // NOLINT
}

template <typename T>
static void declare_operator_atom(nb::module_ &m, std::string const &type_name) {
    std::string pyclass_name = "OperatorAtom" + type_name;
    using basis_t = typename OperatorAtom<T>::basis_t;
    nb::class_<OperatorAtom<T>, Operator<OperatorAtom<T>>> pyclass(m, pyclass_name.c_str());
    pyclass.def(nb::init<std::shared_ptr<const basis_t>>())
        .def(nb::init<std::shared_ptr<const basis_t>, OperatorType, int>());
}

template <typename T>
static void declare_operator_pair(nb::module_ &m, std::string const &type_name) {
    std::string pyclass_name = "OperatorPair" + type_name;
    using basis_t = typename OperatorPair<T>::basis_t;
    nb::class_<OperatorPair<T>, Operator<OperatorPair<T>>> pyclass(m, pyclass_name.c_str());
    pyclass.def(nb::init<std::shared_ptr<const basis_t>>())
        .def(nb::init<std::shared_ptr<const basis_t>, OperatorType>());
}

void bind_operator(nb::module_ &m) {
    declare_operator<OperatorAtom<double>>(m, "OperatorAtomReal");
    declare_operator<OperatorAtom<std::complex<double>>>(m, "OperatorAtomComplex");
    declare_operator_atom<double>(m, "Real");
    declare_operator_atom<std::complex<double>>(m, "Complex");

    declare_operator<OperatorPair<double>>(m, "OperatorPairReal");
    declare_operator<OperatorPair<std::complex<double>>>(m, "OperatorPairComplex");
    declare_operator_pair<double>(m, "Real");
    declare_operator_pair<std::complex<double>>(m, "Complex");
}
