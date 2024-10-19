#include "./Basis.py.hpp"

#include "pairinteraction/basis/Basis.hpp"
#include "pairinteraction/basis/BasisAtom.hpp"
#include "pairinteraction/basis/BasisAtomCreator.hpp"
#include "pairinteraction/basis/BasisClassicalLight.hpp"
#include "pairinteraction/basis/BasisClassicalLightCreator.hpp"
#include "pairinteraction/basis/BasisCombined.hpp"
#include "pairinteraction/basis/BasisCombinedCreator.hpp"
#include "pairinteraction/database/Database.hpp"
#include "pairinteraction/interfaces/TransformationBuilderInterface.hpp"
#include "pairinteraction/ket/KetAtom.hpp"
#include "pairinteraction/ket/KetCombined.hpp"
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
        .def("transformed", nb::overload_cast<const Sorting &>(&Basis<T>::transformed, nb::const_));
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
static void declare_basis_classical_light(nb::module_ &m, std::string const &type_name) {
    std::string pyclass_name = "BasisClassicalLight" + type_name;
    nb::class_<BasisClassicalLight<T>, Basis<BasisClassicalLight<T>>> pyclass(m,
                                                                              pyclass_name.c_str());
}

template <typename T>
static void declare_basis_classical_light_creator(nb::module_ &m, std::string const &type_name) {
    std::string pyclass_name = "BasisClassicalLightCreator" + type_name;
    nb::class_<BasisClassicalLightCreator<T>> pyclass(m, pyclass_name.c_str());
    pyclass.def(nb::init<>())
        .def("set_photon_energy", &BasisClassicalLightCreator<T>::set_photon_energy)
        .def("restrict_quantum_number_q", &BasisClassicalLightCreator<T>::restrict_quantum_number_q)
        .def("create", &BasisClassicalLightCreator<T>::create);
}

template <typename T>
static void declare_basis_combined(nb::module_ &m, std::string const &type_name) {
    std::string pyclass_name = "BasisCombined" + type_name;
    nb::class_<BasisCombined<T>, Basis<BasisCombined<T>>> pyclass(m, pyclass_name.c_str());
}

template <typename T>
static void declare_basis_combined_creator(nb::module_ &m, std::string const &type_name) {
    std::string pyclass_name = "BasisCombinedCreator" + type_name;
    nb::class_<BasisCombinedCreator<T>> pyclass(m, pyclass_name.c_str());
    pyclass.def(nb::init<>())
        .def("add", &BasisCombinedCreator<T>::add)
        .def("restrict_energy", &BasisCombinedCreator<T>::restrict_energy)
        .def("restrict_quantum_number_m", &BasisCombinedCreator<T>::restrict_quantum_number_m)
        .def("create", &BasisCombinedCreator<T>::create);
}

void bind_basis(nb::module_ &m) {
    declare_basis<BasisAtom<float>>(m, "BasisAtomFloat");
    declare_basis<BasisAtom<double>>(m, "BasisAtomDouble");
    declare_basis<BasisAtom<std::complex<float>>>(m, "BasisAtomComplexFloat");
    declare_basis<BasisAtom<std::complex<double>>>(m, "BasisAtomComplexDouble");
    declare_basis_atom<float>(m, "Float");
    declare_basis_atom<double>(m, "Double");
    declare_basis_atom<std::complex<float>>(m, "ComplexFloat");
    declare_basis_atom<std::complex<double>>(m, "ComplexDouble");
    declare_basis_atom_creator<float>(m, "Float");
    declare_basis_atom_creator<double>(m, "Double");
    declare_basis_atom_creator<std::complex<float>>(m, "ComplexFloat");
    declare_basis_atom_creator<std::complex<double>>(m, "ComplexDouble");

    declare_basis<BasisClassicalLight<float>>(m, "BasisClassicalLightFloat");
    declare_basis<BasisClassicalLight<double>>(m, "BasisClassicalLightDouble");
    declare_basis<BasisClassicalLight<std::complex<float>>>(m, "BasisClassicalLightComplexFloat");
    declare_basis<BasisClassicalLight<std::complex<double>>>(m, "BasisClassicalLightComplexDouble");
    declare_basis_classical_light<float>(m, "Float");
    declare_basis_classical_light<double>(m, "Double");
    declare_basis_classical_light<std::complex<float>>(m, "ComplexFloat");
    declare_basis_classical_light<std::complex<double>>(m, "ComplexDouble");
    declare_basis_classical_light_creator<float>(m, "Float");
    declare_basis_classical_light_creator<double>(m, "Double");
    declare_basis_classical_light_creator<std::complex<float>>(m, "ComplexFloat");
    declare_basis_classical_light_creator<std::complex<double>>(m, "ComplexDouble");

    declare_basis<BasisCombined<float>>(m, "BasisCombinedFloat");
    declare_basis<BasisCombined<double>>(m, "BasisCombinedDouble");
    declare_basis<BasisCombined<std::complex<float>>>(m, "BasisCombinedComplexFloat");
    declare_basis<BasisCombined<std::complex<double>>>(m, "BasisCombinedComplexDouble");
    declare_basis_combined<float>(m, "Float");
    declare_basis_combined<double>(m, "Double");
    declare_basis_combined<std::complex<float>>(m, "ComplexFloat");
    declare_basis_combined<std::complex<double>>(m, "ComplexDouble");
    declare_basis_combined_creator<float>(m, "Float");
    declare_basis_combined_creator<double>(m, "Double");
    declare_basis_combined_creator<std::complex<float>>(m, "ComplexFloat");
    declare_basis_combined_creator<std::complex<double>>(m, "ComplexDouble");
}
