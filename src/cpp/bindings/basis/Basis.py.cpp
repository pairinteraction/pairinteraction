// SPDX-FileCopyrightText: 2024 Pairinteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "./Basis.py.hpp"

#include "pairinteraction/basis/Basis.hpp"
#include "pairinteraction/basis/BasisAtom.hpp"
#include "pairinteraction/basis/BasisAtomCreator.hpp"
#include "pairinteraction/basis/BasisPair.hpp"
#include "pairinteraction/basis/BasisPairCreator.hpp"
#include "pairinteraction/database/Database.hpp"
#include "pairinteraction/interfaces/TransformationBuilderInterface.hpp"
#include "pairinteraction/ket/KetAtom.hpp"
#include "pairinteraction/ket/KetPair.hpp"
#include "pairinteraction/system/SystemAtom.hpp"

#include <nanobind/eigen/sparse.h>
#include <nanobind/nanobind.h>
#include <nanobind/stl/complex.h>
#include <nanobind/stl/set.h>
#include <nanobind/stl/shared_ptr.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

namespace nb = nanobind;
using namespace pairinteraction;

template <typename T>
static void declare_basis(nb::module_ &m, std::string const &type_name) {
    std::string pyclass_name = "Basis" + type_name;
    using scalar_t = typename Basis<T>::scalar_t;
    nb::class_<Basis<T>, TransformationBuilderInterface<scalar_t>> pyclass(m, pyclass_name.c_str());
    pyclass.def("get_kets", &Basis<T>::get_kets)
        .def("get_ket", &Basis<T>::get_ket)
        .def("get_state", &Basis<T>::get_state)
        .def("get_number_of_states", &Basis<T>::get_number_of_states)
        .def("get_number_of_kets", &Basis<T>::get_number_of_kets)
        .def("get_quantum_number_f", &Basis<T>::get_quantum_number_f)
        .def("get_quantum_number_m", &Basis<T>::get_quantum_number_m)
        .def("get_parity", &Basis<T>::get_parity)
        .def("get_coefficients", nb::overload_cast<>(&Basis<T>::get_coefficients, nb::const_))
        .def("get_transformation", &Basis<T>::get_transformation)
        .def("get_rotator", &Basis<T>::get_rotator)
        .def("get_sorter", &Basis<T>::get_sorter)
        .def("get_indices_of_blocks", &Basis<T>::get_indices_of_blocks)
        .def("get_sorter_without_checks", &Basis<T>::get_sorter_without_checks)
        .def("get_indices_of_blocks_without_checks",
             &Basis<T>::get_indices_of_blocks_without_checks)
        .def(
            "transformed",
            nb::overload_cast<const Transformation<scalar_t> &>(&Basis<T>::transformed, nb::const_))
        .def("transformed", nb::overload_cast<const Sorting &>(&Basis<T>::transformed, nb::const_))
        .def("get_amplitudes",
             nb::overload_cast<std::shared_ptr<const typename Basis<T>::ket_t>>(
                 &Basis<T>::get_amplitudes, nb::const_))
        .def("get_amplitudes",
             nb::overload_cast<std::shared_ptr<const T>>(&Basis<T>::get_amplitudes, nb::const_))
        .def("get_overlaps",
             nb::overload_cast<std::shared_ptr<const typename Basis<T>::ket_t>>(
                 &Basis<T>::get_overlaps, nb::const_))
        .def("get_overlaps",
             nb::overload_cast<std::shared_ptr<const T>>(&Basis<T>::get_overlaps, nb::const_))
        .def("get_matrix_elements",
             nb::overload_cast<std::shared_ptr<const typename Basis<T>::ket_t>, OperatorType, int>(
                 &Basis<T>::get_matrix_elements, nb::const_))
        .def("get_matrix_elements",
             nb::overload_cast<std::shared_ptr<const T>, OperatorType, int>(
                 &Basis<T>::get_matrix_elements, nb::const_))
        .def("get_corresponding_state",
             nb::overload_cast<size_t>(&Basis<T>::get_corresponding_state, nb::const_))
        .def("get_corresponding_state",
             nb::overload_cast<std::shared_ptr<const typename Basis<T>::ket_t>>(
                 &Basis<T>::get_corresponding_state, nb::const_))
        .def("get_corresponding_state_index",
             nb::overload_cast<size_t>(&Basis<T>::get_corresponding_state_index, nb::const_))
        .def("get_corresponding_state_index",
             nb::overload_cast<std::shared_ptr<const typename Basis<T>::ket_t>>(
                 &Basis<T>::get_corresponding_state_index, nb::const_))
        .def("get_corresponding_ket",
             nb::overload_cast<size_t>(&Basis<T>::get_corresponding_ket, nb::const_))
        .def("get_corresponding_ket",
             nb::overload_cast<std::shared_ptr<const T>>(&Basis<T>::get_corresponding_ket,
                                                         nb::const_))
        .def("get_corresponding_ket_index",
             nb::overload_cast<size_t>(&Basis<T>::get_corresponding_ket_index, nb::const_))
        .def("get_corresponding_ket_index",
             nb::overload_cast<std::shared_ptr<const T>>(&Basis<T>::get_corresponding_ket_index,
                                                         nb::const_));
}

template <typename T>
static void declare_basis_atom(nb::module_ &m, std::string const &type_name) {
    std::string pyclass_name = "BasisAtom" + type_name;
    nb::class_<BasisAtom<T>, Basis<BasisAtom<T>>> pyclass(m, pyclass_name.c_str());
}

template <typename T>
static void declare_basis_atom_creator(nb::module_ &m, std::string const &type_name) {
    std::string pyclass_name = "BasisAtomCreator" + type_name;
    nb::class_<BasisAtomCreator<T>> pyclass(m, pyclass_name.c_str());
    pyclass.def(nb::init<>())
        .def("set_species", &BasisAtomCreator<T>::set_species)
        .def("restrict_energy", &BasisAtomCreator<T>::restrict_energy)
        .def("restrict_quantum_number_f", &BasisAtomCreator<T>::restrict_quantum_number_f)
        .def("restrict_quantum_number_m", &BasisAtomCreator<T>::restrict_quantum_number_m)
        .def("restrict_parity", &BasisAtomCreator<T>::restrict_parity)
        .def("restrict_quantum_number_n", &BasisAtomCreator<T>::restrict_quantum_number_n)
        .def("restrict_quantum_number_nu", &BasisAtomCreator<T>::restrict_quantum_number_nu)
        .def("restrict_quantum_number_l", &BasisAtomCreator<T>::restrict_quantum_number_l)
        .def("restrict_quantum_number_s", &BasisAtomCreator<T>::restrict_quantum_number_s)
        .def("restrict_quantum_number_j", &BasisAtomCreator<T>::restrict_quantum_number_j)
        .def("append_ket", &BasisAtomCreator<T>::append_ket)
        .def("create", &BasisAtomCreator<T>::create);
}

template <typename T>
static void declare_basis_pair(nb::module_ &m, std::string const &type_name) {
    std::string pyclass_name = "BasisPair" + type_name;
    nb::class_<BasisPair<T>, Basis<BasisPair<T>>> pyclass(m, pyclass_name.c_str());
    pyclass
        .def("get_amplitudes",
             nb::overload_cast<std::shared_ptr<const KetAtom>, std::shared_ptr<const KetAtom>>(
                 &BasisPair<T>::get_amplitudes, nb::const_))
        .def("get_amplitudes",
             nb::overload_cast<std::shared_ptr<const BasisAtom<T>>,
                               std::shared_ptr<const BasisAtom<T>>>(&BasisPair<T>::get_amplitudes,
                                                                    nb::const_))
        .def("get_overlaps",
             nb::overload_cast<std::shared_ptr<const KetAtom>, std::shared_ptr<const KetAtom>>(
                 &BasisPair<T>::get_overlaps, nb::const_))
        .def("get_overlaps",
             nb::overload_cast<std::shared_ptr<const BasisAtom<T>>,
                               std::shared_ptr<const BasisAtom<T>>>(&BasisPair<T>::get_overlaps,
                                                                    nb::const_))
        .def("get_matrix_elements",
             nb::overload_cast<std::shared_ptr<const BasisPair<T>>, OperatorType, OperatorType, int,
                               int>(&BasisPair<T>::get_matrix_elements, nb::const_))
        .def("get_matrix_elements",
             nb::overload_cast<std::shared_ptr<const BasisAtom<T>>,
                               std::shared_ptr<const BasisAtom<T>>, OperatorType, OperatorType, int,
                               int>(&BasisPair<T>::get_matrix_elements, nb::const_))
        .def("get_matrix_elements",
             nb::overload_cast<std::shared_ptr<const typename BasisPair<T>::ket_t>, OperatorType,
                               OperatorType, int, int>(&BasisPair<T>::get_matrix_elements,
                                                       nb::const_))
        .def("get_matrix_elements",
             nb::overload_cast<std::shared_ptr<const KetAtom>, std::shared_ptr<const KetAtom>,
                               OperatorType, OperatorType, int, int>(
                 &BasisPair<T>::get_matrix_elements, nb::const_));
}

template <typename T>
static void declare_basis_pair_creator(nb::module_ &m, std::string const &type_name) {
    std::string pyclass_name = "BasisPairCreator" + type_name;
    nb::class_<BasisPairCreator<T>> pyclass(m, pyclass_name.c_str());
    pyclass.def(nb::init<>())
        .def("add", &BasisPairCreator<T>::add)
        .def("restrict_energy", &BasisPairCreator<T>::restrict_energy)
        .def("restrict_quantum_number_m", &BasisPairCreator<T>::restrict_quantum_number_m)
        .def("restrict_product_of_parities", &BasisPairCreator<T>::restrict_product_of_parities)
        .def("create", &BasisPairCreator<T>::create);
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
