#include "./matrix_elements.py.hpp"

#include "pairinteraction/convenience/matrix_elements.hpp"
#include "pairinteraction/ket/KetAtom.hpp"
#include "pairinteraction/system/SystemAtom.hpp"

#include <complex>
#include <nanobind/nanobind.h>

namespace nb = nanobind;
using namespace pairinteraction;

template <typename T>
static void declare_calculate_energy(nb::module_ &m) {
    using ket_ptr_t = std::shared_ptr<const KetAtom<T>>;

    m.def("calculate_energy",
          nb::overload_cast<ket_ptr_t, const SystemAtom<T> &>(&calculate_energy<T>));
    m.def("calculate_energy",
          nb::overload_cast<ket_ptr_t, const SystemAtom<std::complex<T>> &>(
              &calculate_energy<std::complex<T>>));
    m.def("calculate_energy", nb::overload_cast<ket_ptr_t>(&calculate_energy<T>));
}

template <typename T>
static void declare_calculate_electric_dipole_matrix_element(nb::module_ &m) {
    using ket_ptr_t = std::shared_ptr<const KetAtom<T>>;

    m.def("declare_calculate_electric_dipole_matrix_element",
          nb::overload_cast<ket_ptr_t, ket_ptr_t, const SystemAtom<T> &, int>(
              &calculate_electric_dipole_matrix_element<T>));
    m.def("declare_calculate_electric_dipole_matrix_element",
          nb::overload_cast<ket_ptr_t, ket_ptr_t, const SystemAtom<std::complex<T>> &, int>(
              &calculate_electric_dipole_matrix_element<std::complex<T>>));
    m.def(
        "declare_calculate_electric_dipole_matrix_element",
        nb::overload_cast<ket_ptr_t, ket_ptr_t, int>(&calculate_electric_dipole_matrix_element<T>));
}

void bind_matrix_elements(nb::module_ &m) {
    declare_calculate_energy<float>(m);
    declare_calculate_energy<double>(m);
    declare_calculate_electric_dipole_matrix_element<float>(m);
    declare_calculate_electric_dipole_matrix_element<double>(m);
}
