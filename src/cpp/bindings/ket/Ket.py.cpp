#include "./Ket.py.hpp"

#include "pairinteraction/database/Database.hpp"
#include "pairinteraction/ket/Ket.hpp"
#include "pairinteraction/ket/KetAtom.hpp"
#include "pairinteraction/ket/KetAtomCreator.hpp"
#include "pairinteraction/ket/KetClassicalLight.hpp"
#include "pairinteraction/ket/KetClassicalLightCreator.hpp"
#include "pairinteraction/ket/KetCombined.hpp"
#include "pairinteraction/utils/traits.hpp"

#include <complex>
#include <nanobind/nanobind.h>
#include <nanobind/stl/shared_ptr.h>
#include <nanobind/stl/string.h>
#include <sstream>

namespace nb = nanobind;
using namespace pairinteraction;

template <typename T>
static void declare_ket(nb::module_ &m, std::string const &type_name) {
    std::string pyclass_name = "Ket" + type_name;
    nb::class_<Ket<T>> pyclass(m, pyclass_name.c_str());
    pyclass.def("get_energy", &Ket<T>::get_energy)
        .def("get_quantum_number_f", &Ket<T>::get_quantum_number_f)
        .def("get_quantum_number_m", &Ket<T>::get_quantum_number_m)
        .def("get_parity", &Ket<T>::get_parity)
        .def("get_label", &Ket<T>::get_label)
        .def("__str__", [](Ket<T> const &self) {
            std::stringstream ss;
            ss << self;
            return ss.str();
        });
}

template <typename T>
static void declare_ket_atom(nb::module_ &m, std::string const &type_name) {
    std::string pyclass_name = "KetAtom" + type_name;
    nb::class_<KetAtom<T>, Ket<T>> pyclass(m, pyclass_name.c_str());
    pyclass.def("get_species", &KetAtom<T>::get_species)
        .def("get_quantum_number_n", &KetAtom<T>::get_quantum_number_n)
        .def("get_quantum_number_nu", &KetAtom<T>::get_quantum_number_nu)
        .def("get_quantum_number_l", &KetAtom<T>::get_quantum_number_l)
        .def("get_quantum_number_s", &KetAtom<T>::get_quantum_number_s)
        .def("get_quantum_number_j", &KetAtom<T>::get_quantum_number_j);
}

template <typename T>
static void declare_ket_atom_creator(nb::module_ &m, std::string const &type_name) {
    std::string pyclass_name = "KetAtomCreator" + type_name;
    nb::class_<KetAtomCreator<T>> pyclass(m, pyclass_name.c_str());
    pyclass.def(nb::init<>())
        .def(nb::init<std::string, int, T, T, T>())
        .def("set_species", &KetAtomCreator<T>::set_species)
        .def("set_energy", &KetAtomCreator<T>::set_energy)
        .def("set_quantum_number_f", &KetAtomCreator<T>::set_quantum_number_f)
        .def("set_quantum_number_m", &KetAtomCreator<T>::set_quantum_number_m)
        .def("set_parity", &KetAtomCreator<T>::set_parity)
        .def("set_quantum_number_n", &KetAtomCreator<T>::set_quantum_number_n)
        .def("set_quantum_number_nu", &KetAtomCreator<T>::set_quantum_number_nu)
        .def("set_quantum_number_l", &KetAtomCreator<T>::set_quantum_number_l)
        .def("set_quantum_number_s", &KetAtomCreator<T>::set_quantum_number_s)
        .def("set_quantum_number_j", &KetAtomCreator<T>::set_quantum_number_j)
        .def("create", &KetAtomCreator<T>::create);
}

template <typename T>
static void declare_ket_classical_light(nb::module_ &m, std::string const &type_name) {
    std::string pyclass_name = "KetClassicalLight" + type_name;
    nb::class_<KetClassicalLight<T>, Ket<T>> pyclass(m, pyclass_name.c_str());
    pyclass.def("get_photon_energy", &KetClassicalLight<T>::get_photon_energy)
        .def("get_quantum_number_q", &KetClassicalLight<T>::get_quantum_number_q);
}

template <typename T>
static void declare_ket_classical_light_creator(nb::module_ &m, std::string const &type_name) {
    std::string pyclass_name = "KetClassicalLightCreator" + type_name;
    nb::class_<KetClassicalLightCreator<T>> pyclass(m, pyclass_name.c_str());
    pyclass.def(nb::init<>())
        .def(nb::init<T, int>())
        .def("set_photon_energy", &KetClassicalLightCreator<T>::set_photon_energy)
        .def("set_quantum_number_q", &KetClassicalLightCreator<T>::set_quantum_number_q)
        .def("create", &KetClassicalLightCreator<T>::create);
}

template <typename T>
static void declare_ket_combined(nb::module_ &m, std::string const &type_name) {
    std::string pyclass_name = "KetCombined" + type_name;
    using real_t = typename traits::NumTraits<T>::real_t;
    nb::class_<KetCombined<T>, Ket<real_t>> pyclass(m, pyclass_name.c_str());
}

void bind_ket(nb::module_ &m) {
    declare_ket<float>(m, "Float");
    declare_ket<double>(m, "Double");
    declare_ket_atom<float>(m, "Float");
    declare_ket_atom<double>(m, "Double");
    declare_ket_atom_creator<float>(m, "Float");
    declare_ket_atom_creator<double>(m, "Double");
    declare_ket_classical_light<float>(m, "Float");
    declare_ket_classical_light<double>(m, "Double");
    declare_ket_classical_light_creator<float>(m, "Float");
    declare_ket_classical_light_creator<double>(m, "Double");
    declare_ket_combined<float>(m, "Float");
    declare_ket_combined<double>(m, "Double");
    declare_ket_combined<std::complex<float>>(m, "ComplexFloat");
    declare_ket_combined<std::complex<double>>(m, "ComplexDouble");
}
