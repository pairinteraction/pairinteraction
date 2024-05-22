#include "Basis.py.hpp"

#include "Basis.hpp"
#include "BasisAtom.hpp"
#include "BasisAtomCreator.hpp"
#include "BasisClassicalLight.hpp"
#include "BasisClassicalLightCreator.hpp"
#include "database/Database.hpp"
#include "interfaces/TransformationBuilderInterface.hpp"

#include <nanobind/nanobind.h>
#include <nanobind/eigen/sparse.h>
#include <nanobind/stl/complex.h>
#include <nanobind/stl/vector.h>

namespace nb = nanobind;

template <typename T>
static void declare_basis(nb::module_ &m, std::string type_name) {
    std::string pylass_name = "Basis" + type_name;
    using scalar_t = typename Basis<T>::scalar_t;
    nb::class_<Basis<T>, TransformationBuilderInterface<scalar_t>> pyclass(m, pylass_name.c_str());
    pyclass
        .def("get_kets", &Basis<T>::get_kets)
        .def("get_number_of_states", &Basis<T>::get_number_of_states)
        .def("get_number_of_kets", &Basis<T>::get_number_of_kets)
        .def("get_quantum_number_f", &Basis<T>::get_quantum_number_f)
        .def("get_quantum_number_m", &Basis<T>::get_quantum_number_m)
        .def("get_parity", &Basis<T>::get_parity)
        .def("get_coefficients", &Basis<T>::get_coefficients)
        .def("get_transformation", &Basis<T>::get_transformation)
        .def("get_rotator", &Basis<T>::get_rotator)
        .def("get_sorter", &Basis<T>::get_sorter)
        .def("get_blocks", &Basis<T>::get_blocks)
        .def("get_sorter_without_checks", &Basis<T>::get_sorter_without_checks)
        .def("get_blocks_without_checks", &Basis<T>::get_blocks_without_checks)
        .def("transform", nb::overload_cast<const Transformation<scalar_t> &>(&Basis<T>::transform, nb::const_))
        .def("transform", nb::overload_cast<const Sorting &>(&Basis<T>::transform, nb::const_));
}

template <typename T>
static void declare_basis_atom(nb::module_ &m, std::string type_name) {
    std::string pylass_name = "BasisAtom" + type_name;
    nb::class_<BasisAtom<T>, Basis<BasisAtom<T>>> pyclass(m, pylass_name.c_str());
    pyclass
        .def("get_database", &BasisAtom<T>::get_database);
}

template <typename T>
static void declare_basis_atom_creator(nb::module_ &m, std::string type_name) {
    std::string pylass_name = "BasisAtomCreator" + type_name;
    nb::class_<BasisAtomCreator<T>> pyclass(m, pylass_name.c_str());
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
        .def("add_ket", &BasisAtomCreator<T>::add_ket)
        .def("create", [](BasisAtomCreator<T> const &basis_atom_creator) {
            Database &database = Database::get_global_instance();
            basis_atom_creator.create(database);
        });
}

template <typename T>
static void declare_basis_classical_light(nb::module_ &m, std::string type_name) {
    std::string pylass_name = "BasisClassicalLight" + type_name;
    nb::class_<BasisClassicalLight<T>, Basis<BasisClassicalLight<T>>> pyclass(m, pylass_name.c_str());
}

template <typename T>
static void declare_basis_classical_light_creator(nb::module_ &m, std::string type_name) {
    std::string pylass_name = "BasisClassicalLightCreator" + type_name;
    nb::class_<BasisClassicalLightCreator<T>> pyclass(m, pylass_name.c_str());
    pyclass.def(nb::init<>())
        .def("set_photon_energy", &BasisClassicalLightCreator<T>::set_photon_energy)
        .def("restrict_quantum_number_q", &BasisClassicalLightCreator<T>::restrict_quantum_number_q)
        .def("create", &BasisClassicalLightCreator<T>::create);
}

void bind_basis(nb::module_ &m) {
    declare_basis<BasisAtom<float>>(m, "BasisAtomFloat");
    declare_basis<BasisAtom<double>>(m, "BasisAtomDouble");
    //declare_basis<BasisAtom<std::complex<float>>>(m, "BasisAtomComplexFloat");
    //declare_basis<BasisAtom<std::complex<double>>>(m, "BasisAtomComplexDouble");
    declare_basis<BasisClassicalLight<float>>(m, "BasisClassicalLightFloat");
    declare_basis<BasisClassicalLight<double>>(m, "BasisClassicalLightDouble");
    //declare_basis<BasisClassicalLight<std::complex<float>>>(m, "BasisClassicalLightComplexFloat");
    //declare_basis<BasisClassicalLight<std::complex<double>>>(m, "BasisClassicalLightComplexDouble");
    declare_basis_atom<float>(m, "Float");
    declare_basis_atom<double>(m, "Double");
    //declare_basis_atom<std::complex<float>>(m, "ComplexFloat");
    //declare_basis_atom<std::complex<double>>(m, "ComplexDouble");
    declare_basis_atom_creator<float>(m, "Float");
    declare_basis_atom_creator<double>(m, "Double");
    //declare_basis_atom_creator<std::complex<float>>(m, "ComplexFloat");
    //declare_basis_atom_creator<std::complex<double>>(m, "ComplexDouble");
    declare_basis_classical_light<float>(m, "Float");
    declare_basis_classical_light<double>(m, "Double");
    //declare_basis_classical_light<std::complex<float>>(m, "ComplexFloat");
    //declare_basis_classical_light<std::complex<double>>(m, "ComplexDouble");
    declare_basis_classical_light_creator<float>(m, "Float");
    declare_basis_classical_light_creator<double>(m, "Double");
    //declare_basis_classical_light_creator<std::complex<float>>(m, "ComplexFloat");
    //declare_basis_classical_light_creator<std::complex<double>>(m, "ComplexDouble");
}
