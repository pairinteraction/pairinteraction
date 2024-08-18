#include "Ket.py.hpp"

#include "pairinteraction/database/Database.hpp"
#include "pairinteraction/ket/Ket.hpp"
#include "pairinteraction/ket/KetAtom.hpp"
#include "pairinteraction/ket/KetAtomCreator.hpp"

#include <nanobind/nanobind.h>
#include <sstream>

namespace nb = nanobind;

template <typename T>
static void declare_ket(nb::module_ &m, const std::string &type_name) {
    std::string pylass_name = "Ket" + type_name;
    nb::class_<pairinteraction::Ket<T>> pyclass(m, pylass_name.c_str());
    pyclass.def("get_energy", &pairinteraction::Ket<T>::get_energy)
        .def("get_quantum_number_f", &pairinteraction::Ket<T>::get_quantum_number_f)
        .def("get_quantum_number_m", &pairinteraction::Ket<T>::get_quantum_number_m)
        .def("get_parity", &pairinteraction::Ket<T>::get_parity)
        .def("get_label", &pairinteraction::Ket<T>::get_label)
        .def("__str__", [](pairinteraction::Ket<T> const &self) {
            std::stringstream ss;
            ss << self;
            return ss.str();
        });
}

template <typename T>
static void declare_ket_atom(nb::module_ &m, const std::string &type_name) {
    std::string pylass_name = "KetAtom" + type_name;
    nb::class_<pairinteraction::KetAtom<T>, pairinteraction::Ket<T>> pyclass(m,
                                                                             pylass_name.c_str());
    pyclass.def("get_species", &pairinteraction::KetAtom<T>::get_species)
        .def("get_quantum_number_n", &pairinteraction::KetAtom<T>::get_quantum_number_n)
        .def("get_quantum_number_nu", &pairinteraction::KetAtom<T>::get_quantum_number_nu)
        .def("get_quantum_number_l", &pairinteraction::KetAtom<T>::get_quantum_number_l)
        .def("get_quantum_number_s", &pairinteraction::KetAtom<T>::get_quantum_number_s)
        .def("get_quantum_number_j", &pairinteraction::KetAtom<T>::get_quantum_number_j);
}

template <typename T>
static void declare_ket_atom_creator(nb::module_ &m, const std::string &type_name) {
    std::string pylass_name = "KetAtomCreator" + type_name;
    nb::class_<pairinteraction::KetAtomCreator<T>> pyclass(m, pylass_name.c_str());
    pyclass.def(nb::init<>())
        .def(nb::init<std::string, int, T, T, T>())
        .def("set_species", &pairinteraction::KetAtomCreator<T>::set_species)
        .def("set_energy", &pairinteraction::KetAtomCreator<T>::set_energy)
        .def("set_quantum_number_f", &pairinteraction::KetAtomCreator<T>::set_quantum_number_f)
        .def("set_quantum_number_m", &pairinteraction::KetAtomCreator<T>::set_quantum_number_m)
        .def("set_parity", &pairinteraction::KetAtomCreator<T>::set_parity)
        .def("set_quantum_number_n", &pairinteraction::KetAtomCreator<T>::set_quantum_number_n)
        .def("set_quantum_number_nu", &pairinteraction::KetAtomCreator<T>::set_quantum_number_nu)
        .def("set_quantum_number_l", &pairinteraction::KetAtomCreator<T>::set_quantum_number_l)
        .def("set_quantum_number_s", &pairinteraction::KetAtomCreator<T>::set_quantum_number_s)
        .def("set_quantum_number_j", &pairinteraction::KetAtomCreator<T>::set_quantum_number_j)
        .def("create", [](pairinteraction::KetAtomCreator<T> const &self) {
            thread_local static pairinteraction::Database db;
            return self.create(db);
        });
}

void bind_ket(nb::module_ &m) {
    declare_ket<float>(m, "Float");
    declare_ket<double>(m, "Double");
    declare_ket_atom<float>(m, "Float");
    declare_ket_atom<double>(m, "Double");
    declare_ket_atom_creator<float>(m, "Float");
    declare_ket_atom_creator<double>(m, "Double");
}
