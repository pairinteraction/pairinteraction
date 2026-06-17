// SPDX-FileCopyrightText: 2024 PairInteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "./Ket.py.hpp"

#include "pairinteraction/basis/BasisAtom.hpp"
#include "pairinteraction/database/Database.hpp"
#include "pairinteraction/ket/Ket.hpp"
#include "pairinteraction/ket/KetAtom.hpp"
#include "pairinteraction/ket/KetAtomCreator.hpp"
#include "pairinteraction/ket/KetPair.hpp"
#include "pairinteraction/utils/traits.hpp"

#include <complex>
#include <nanobind/nanobind.h>
#include <nanobind/operators.h>
#include <nanobind/stl/shared_ptr.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>
#include <sstream>

namespace nb = nanobind;
using namespace pairinteraction;

static void declare_ket(nb::module_ &m) {
    std::string pyclass_name = "Ket";
    nb::class_<Ket> pyclass(m, pyclass_name.c_str());
    pyclass.def("get_energy", &Ket::get_energy)
        .def("get_label", &Ket::get_label)
        .def("__str__", [](Ket const &self) {
            std::stringstream ss;
            ss << self;
            return ss.str();
        });
}

static void declare_ket_atom(nb::module_ &m) {
    std::string pyclass_name = "KetAtom";
    nb::class_<KetAtom, Ket> pyclass(m, pyclass_name.c_str());
    pyclass.def("get_database", &KetAtom::get_database, nb::rv_policy::reference)
        .def("get_species", &KetAtom::get_species)
        .def("get_quantum_number", &KetAtom::get_quantum_number)
        .def("get_quantum_number_std", &KetAtom::get_quantum_number_std)
        .def(nb::self == nb::self) // NOLINT(misc-redundant-expression)
        .def("__hash__", [](const KetAtom &self) { return KetAtom::hash{}(self); });
}

static void declare_ket_atom_creator(nb::module_ &m) {
    std::string pyclass_name = "KetAtomCreator";
    nb::class_<KetAtomCreator> pyclass(m, pyclass_name.c_str());
    pyclass.def(nb::init<>())
        .def(nb::init<std::string, int, double, double, double>())
        .def("set_species", &KetAtomCreator::set_species)
        .def("set_energy", &KetAtomCreator::set_energy)
        .def("set_quantum_number", &KetAtomCreator::set_quantum_number)
        .def("create", &KetAtomCreator::create);
}

template <typename T>
static void declare_ket_pair(nb::module_ &m, std::string const &type_name) {
    std::string pyclass_name = "KetPair" + type_name;
    nb::class_<KetPair<T>, Ket> pyclass(m, pyclass_name.c_str());
    pyclass.def("get_atomic_states", &KetPair<T>::get_atomic_states)
        .def("get_quantum_number_m", &KetPair<T>::get_quantum_number_m)
        .def(nb::self == nb::self) // NOLINT(misc-redundant-expression)
        .def("__hash__", [](const KetPair<T> &self) { return typename KetPair<T>::hash{}(self); });
}

void bind_ket(nb::module_ &m) {
    declare_ket(m);
    declare_ket_atom(m);
    declare_ket_atom_creator(m);
    declare_ket_pair<double>(m, "Real");
    declare_ket_pair<std::complex<double>>(m, "Complex");
}
