#include "Ket.py.hpp"

#include "Ket.hpp"
#include "KetAtom.hpp"
#include "KetAtomCreator.hpp"

#include <sstream>

#include <nanobind/nanobind.h>

namespace nb = nanobind;

template <typename T>
static void declare_ket(nb::module_ &m, std::string type_name) {
    std::string pylass_name = "Ket" + type_name;
    nb::class_<Ket<T>> pyclass(m, pylass_name.c_str());
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
static void declare_ket_atom(nb::module_ &m, std::string type_name) {
    std::string pylass_name = "KetAtom" + type_name;
    nb::class_<KetAtom<T>, Ket<T>> pyclass(m, pylass_name.c_str());
    pyclass.def("get_species", &KetAtom<T>::get_species)
        .def("get_quantum_number_n", &KetAtom<T>::get_quantum_number_n)
        .def("get_quantum_number_nu", &KetAtom<T>::get_quantum_number_nu)
        .def("get_quantum_number_l", &KetAtom<T>::get_quantum_number_l)
        .def("get_quantum_number_s", &KetAtom<T>::get_quantum_number_s)
        .def("get_quantum_number_j", &KetAtom<T>::get_quantum_number_j);
}

template <typename T>
static void declare_ket_atom_creator(nb::module_ &m, std::string type_name) {
    std::string pylass_name = "KetAtomCreator" + type_name;
    nb::class_<KetAtomCreator<T>> pyclass(m, pylass_name.c_str());
    pyclass.def(nb::init<std::string>())
        .def(nb::init<std::string, int, T, float, float>())
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

void bind_ket(nb::module_ &m) {
    declare_ket<float>(m, "Float");
    declare_ket<double>(m, "Double");
    declare_ket_atom<float>(m, "Float");
    declare_ket_atom<double>(m, "Double");
    declare_ket_atom_creator<float>(m, "Float");
    declare_ket_atom_creator<double>(m, "Double");
}
