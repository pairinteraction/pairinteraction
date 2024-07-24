#include "pintr/ket/Ket.py.hpp"

#include "pintr/database/Database.hpp"
#include "pintr/ket/Ket.hpp"
#include "pintr/ket/KetAtom.hpp"
#include "pintr/ket/KetAtomCreator.hpp"

#include <nanobind/nanobind.h>
#include <sstream>

namespace nb = nanobind;

template <typename T>
static void declare_ket(nb::module_ &m, const std::string &type_name) {
    std::string pylass_name = "Ket" + type_name;
    nb::class_<pintr::Ket<T>> pyclass(m, pylass_name.c_str());
    pyclass.def("get_energy", &pintr::Ket<T>::get_energy)
        .def("get_quantum_number_f", &pintr::Ket<T>::get_quantum_number_f)
        .def("get_quantum_number_m", &pintr::Ket<T>::get_quantum_number_m)
        .def("get_parity", &pintr::Ket<T>::get_parity)
        .def("get_label", &pintr::Ket<T>::get_label)
        .def("__str__", [](pintr::Ket<T> const &self) {
            std::stringstream ss;
            ss << self;
            return ss.str();
        });
}

template <typename T>
static void declare_ket_atom(nb::module_ &m, const std::string &type_name) {
    std::string pylass_name = "KetAtom" + type_name;
    nb::class_<pintr::KetAtom<T>, pintr::Ket<T>> pyclass(m, pylass_name.c_str());
    pyclass.def("get_species", &pintr::KetAtom<T>::get_species)
        .def("get_quantum_number_n", &pintr::KetAtom<T>::get_quantum_number_n)
        .def("get_quantum_number_nu", &pintr::KetAtom<T>::get_quantum_number_nu)
        .def("get_quantum_number_l", &pintr::KetAtom<T>::get_quantum_number_l)
        .def("get_quantum_number_s", &pintr::KetAtom<T>::get_quantum_number_s)
        .def("get_quantum_number_j", &pintr::KetAtom<T>::get_quantum_number_j);
}

template <typename T>
static void declare_ket_atom_creator(nb::module_ &m, const std::string &type_name) {
    std::string pylass_name = "KetAtomCreator" + type_name;
    nb::class_<pintr::KetAtomCreator<T>> pyclass(m, pylass_name.c_str());
    pyclass.def(nb::init<>())
        .def(nb::init<std::string, int, T, T, T>())
        .def("set_species", &pintr::KetAtomCreator<T>::set_species)
        .def("set_energy", &pintr::KetAtomCreator<T>::set_energy)
        .def("set_quantum_number_f", &pintr::KetAtomCreator<T>::set_quantum_number_f)
        .def("set_quantum_number_m", &pintr::KetAtomCreator<T>::set_quantum_number_m)
        .def("set_parity", &pintr::KetAtomCreator<T>::set_parity)
        .def("set_quantum_number_n", &pintr::KetAtomCreator<T>::set_quantum_number_n)
        .def("set_quantum_number_nu", &pintr::KetAtomCreator<T>::set_quantum_number_nu)
        .def("set_quantum_number_l", &pintr::KetAtomCreator<T>::set_quantum_number_l)
        .def("set_quantum_number_s", &pintr::KetAtomCreator<T>::set_quantum_number_s)
        .def("set_quantum_number_j", &pintr::KetAtomCreator<T>::set_quantum_number_j)
        .def("create", [](pintr::KetAtomCreator<T> const &self) {
            thread_local static pintr::Database db;
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
