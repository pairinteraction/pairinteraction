// SPDX-FileCopyrightText: 2024 PairInteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "./Basis.py.hpp"

#include "pairinteraction/basis/Basis.hpp"
#include "pairinteraction/basis/BasisAtom.hpp"
#include "pairinteraction/basis/BasisAtomCreator.hpp"
#include "pairinteraction/basis/BasisPair.hpp"
#include "pairinteraction/basis/BasisPairCreator.hpp"
#include "pairinteraction/database/Database.hpp"
#include "pairinteraction/ket/KetAtom.hpp"
#include "pairinteraction/ket/KetPair.hpp"

#include <nanobind/eigen/sparse.h>
#include <nanobind/nanobind.h>
#include <nanobind/stl/complex.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/shared_ptr.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/unordered_map.h>
#include <nanobind/stl/vector.h>
#include <vector>

namespace nb = nanobind;
using namespace pairinteraction;

template <typename T>
static void declare_basis(nb::module_ &m, std::string const &type_name) {
    std::string pyclass_name = "Basis" + type_name;
    nb::class_<Basis<T>> pyclass(m, pyclass_name.c_str());
    pyclass.def("get_kets", &Basis<T>::get_kets)
        .def("get_ket", &Basis<T>::get_ket)
        .def("get_state", &Basis<T>::get_state)
        .def("get_number_of_states", &Basis<T>::get_number_of_states)
        .def("get_number_of_kets", &Basis<T>::get_number_of_kets)
        .def("has_quantum_number", &Basis<T>::has_quantum_number)
        .def("has_energy", &Basis<T>::has_energy)
        .def("is_sorted_by_energy", &Basis<T>::is_sorted_by_energy)
        .def("get_quantum_number", &Basis<T>::get_quantum_number)
        .def("get_energy", &Basis<T>::get_energy)
        .def("get_coefficients", &Basis<T>::get_coefficients)
        .def("get_modified_copy",
             [](const Basis<T> &self,
                Eigen::SparseMatrix<typename Basis<T>::scalar_t, Eigen::RowMajor> coefficients,
                typename Basis<T>::energies_t energy_of_states,
                typename Basis<T>::quantum_numbers_t quantum_numbers_of_states) {
                 return static_cast<const T &>(self).get_modified_copy(
                     std::move(coefficients), std::move(energy_of_states),
                     std::move(quantum_numbers_of_states));
             })
        .def("get_sorter", &Basis<T>::get_sorter)
        .def("get_sorter_by_energy", &Basis<T>::get_sorter_by_energy)
        .def("get_indices_of_blocks", &Basis<T>::get_indices_of_blocks)
        .def("get_sorted_copy", &Basis<T>::get_sorted_copy)
        .def("get_canonicalized_copy", &Basis<T>::get_canonicalized_copy)
        .def("is_canonical", &Basis<T>::is_canonical)
        .def("merge", [](const Basis<T> &self, const std::shared_ptr<const T> &other) {
            return static_cast<const T &>(self).merge(other);
        });
}

template <typename T>
static void declare_basis_atom(nb::module_ &m, std::string const &type_name) {
    std::string pyclass_name = "BasisAtom" + type_name;
    nb::class_<BasisAtom<T>, Basis<BasisAtom<T>>> pyclass(m, pyclass_name.c_str());
    pyclass.def("get_matrix_elements", &BasisAtom<T>::get_matrix_elements);
}

template <typename T>
static void declare_basis_atom_creator(nb::module_ &m, std::string const &type_name) {
    std::string pyclass_name = "BasisAtomCreator" + type_name;
    nb::class_<BasisAtomCreator<T>> pyclass(m, pyclass_name.c_str());
    // keep_alive because the creator stores a reference to the database
    pyclass.def(nb::init<Database &>(), nb::keep_alive<1, 2>())
        .def("set_species", &BasisAtomCreator<T>::set_species)
        .def("restrict_energy", &BasisAtomCreator<T>::restrict_energy)
        .def("restrict_quantum_number", &BasisAtomCreator<T>::restrict_quantum_number)
        .def("set_quantum_number_standard_deviation_factor",
             &BasisAtomCreator<T>::set_quantum_number_standard_deviation_factor)
        .def("add_ket", &BasisAtomCreator<T>::add_ket)
        .def("create", &BasisAtomCreator<T>::create, nb::call_guard<nb::gil_scoped_release>());
}

template <typename T>
static void declare_basis_pair(nb::module_ &m, std::string const &type_name) {
    std::string pyclass_name = "BasisPair" + type_name;
    nb::class_<BasisPair<T>, Basis<BasisPair<T>>> pyclass(m, pyclass_name.c_str());
    pyclass.def("get_matrix_elements", &BasisPair<T>::get_matrix_elements)
        .def("get_basis1", &BasisPair<T>::get_basis1)
        .def("get_basis2", &BasisPair<T>::get_basis2);
}

template <typename T>
static void declare_basis_pair_creator(nb::module_ &m, std::string const &type_name) {
    std::string pyclass_name = "BasisPairCreator" + type_name;
    nb::class_<BasisPairCreator<T>> pyclass(m, pyclass_name.c_str());
    pyclass
        .def(nb::init<std::shared_ptr<const BasisAtom<T>>, std::shared_ptr<const BasisAtom<T>>>())
        .def("restrict_energy", &BasisPairCreator<T>::restrict_energy)
        .def("restrict_quantum_number", &BasisPairCreator<T>::restrict_quantum_number)
        .def("set_symmetrization_enabled", &BasisPairCreator<T>::set_symmetrization_enabled)
        .def("create", &BasisPairCreator<T>::create, nb::call_guard<nb::gil_scoped_release>());
}

void bind_basis(nb::module_ &m) {
    declare_basis<BasisAtom<double>>(m, "BasisAtomReal");
    declare_basis<BasisAtom<std::complex<double>>>(m, "BasisAtomComplex");
    declare_basis_atom<double>(m, "Real");
    declare_basis_atom<std::complex<double>>(m, "Complex");
    declare_basis_atom_creator<double>(m, "Real");
    declare_basis_atom_creator<std::complex<double>>(m, "Complex");

    declare_basis<BasisPair<double>>(m, "BasisPairReal");
    declare_basis<BasisPair<std::complex<double>>>(m, "BasisPairComplex");
    declare_basis_pair<double>(m, "Real");
    declare_basis_pair<std::complex<double>>(m, "Complex");
    declare_basis_pair_creator<double>(m, "Real");
    declare_basis_pair_creator<std::complex<double>>(m, "Complex");
}
