#include "./Ket.py.hpp"

#include "pairinteraction/basis/BasisAtom.hpp"
#include "pairinteraction/database/Database.hpp"
#include "pairinteraction/ket/Ket.hpp"
#include "pairinteraction/ket/KetAtom.hpp"
#include "pairinteraction/ket/KetAtomCreator.hpp"
#include "pairinteraction/ket/KetClassicalLight.hpp"
#include "pairinteraction/ket/KetClassicalLightCreator.hpp"
#include "pairinteraction/ket/KetPair.hpp"
#include "pairinteraction/utils/traits.hpp"

#include <complex>
#include <nanobind/nanobind.h>
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
        .def("get_quantum_number_f", &Ket::get_quantum_number_f)
        .def("get_quantum_number_m", &Ket::get_quantum_number_m)
        .def("get_parity", &Ket::get_parity)
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
    pyclass.def("get_species", &KetAtom::get_species)
        .def("get_quantum_number_n", &KetAtom::get_quantum_number_n)
        .def("get_quantum_number_nu", &KetAtom::get_quantum_number_nu)
        .def("get_quantum_number_nui", &KetAtom::get_quantum_number_nui)
        .def("get_quantum_number_l", &KetAtom::get_quantum_number_l)
        .def("get_quantum_number_s", &KetAtom::get_quantum_number_s)
        .def("get_quantum_number_j", &KetAtom::get_quantum_number_j)
        .def("get_quantum_number_l_ryd", &KetAtom::get_quantum_number_l_ryd)
        .def("get_quantum_number_j_ryd", &KetAtom::get_quantum_number_j_ryd);
}

static void declare_ket_atom_creator(nb::module_ &m) {
    std::string pyclass_name = "KetAtomCreator";
    nb::class_<KetAtomCreator> pyclass(m, pyclass_name.c_str());
    pyclass.def(nb::init<>())
        .def(nb::init<std::string, int, double, double, double>())
        .def("set_species", &KetAtomCreator::set_species)
        .def("set_energy", &KetAtomCreator::set_energy)
        .def("set_quantum_number_f", &KetAtomCreator::set_quantum_number_f)
        .def("set_quantum_number_m", &KetAtomCreator::set_quantum_number_m)
        .def("set_parity", &KetAtomCreator::set_parity)
        .def("set_quantum_number_n", &KetAtomCreator::set_quantum_number_n)
        .def("set_quantum_number_nu", &KetAtomCreator::set_quantum_number_nu)
        .def("set_quantum_number_nui", &KetAtomCreator::set_quantum_number_nui)
        .def("set_quantum_number_l", &KetAtomCreator::set_quantum_number_l)
        .def("set_quantum_number_s", &KetAtomCreator::set_quantum_number_s)
        .def("set_quantum_number_j", &KetAtomCreator::set_quantum_number_j)
        .def("set_quantum_number_l_ryd", &KetAtomCreator::set_quantum_number_l_ryd)
        .def("set_quantum_number_j_ryd", &KetAtomCreator::set_quantum_number_j_ryd)
        .def("create", &KetAtomCreator::create);
}

static void declare_ket_classical_light(nb::module_ &m) {
    std::string pyclass_name = "KetClassicalLight";
    nb::class_<KetClassicalLight, Ket> pyclass(m, pyclass_name.c_str());
    pyclass.def("get_photon_energy", &KetClassicalLight::get_photon_energy)
        .def("get_quantum_number_q", &KetClassicalLight::get_quantum_number_q);
}

static void declare_ket_classical_light_creator(nb::module_ &m) {
    std::string pyclass_name = "KetClassicalLightCreator";
    nb::class_<KetClassicalLightCreator> pyclass(m, pyclass_name.c_str());
    pyclass.def(nb::init<>())
        .def(nb::init<double, int>())
        .def("set_photon_energy", &KetClassicalLightCreator::set_photon_energy)
        .def("set_quantum_number_q", &KetClassicalLightCreator::set_quantum_number_q)
        .def("create", &KetClassicalLightCreator::create);
}

template <typename T>
static void declare_ket_pair(nb::module_ &m, std::string const &type_name) {
    std::string pyclass_name = "KetPair" + type_name;
    nb::class_<KetPair<T>, Ket> pyclass(m, pyclass_name.c_str());
    pyclass.def("get_atomic_states", &KetPair<T>::get_atomic_states);
}

void bind_ket(nb::module_ &m) {
    declare_ket(m);
    declare_ket_atom(m);
    declare_ket_atom_creator(m);
    declare_ket_classical_light(m);
    declare_ket_classical_light_creator(m);
    declare_ket_pair<double>(m, "Double");
    declare_ket_pair<std::complex<double>>(m, "ComplexDouble");
}
