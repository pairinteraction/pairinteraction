#include "OperatorType.py.hpp"

#include "OperatorType.hpp"

#include <nanobind/nanobind.h>

namespace nb = nanobind;

void bind_operator_type(nb::module_ &m) {
    nb::enum_<OperatorType>(m, "OperatorType")
        .value("ZERO", OperatorType::ZERO)
        .value("ENERGY", OperatorType::ENERGY)
        .value("ELECTRIC_DIPOLE", OperatorType::ELECTRIC_DIPOLE)
        .value("ELECTRIC_QUADRUPOLE", OperatorType::ELECTRIC_QUADRUPOLE)
        .value("ELECTRIC_OCTUPOLE", OperatorType::ELECTRIC_OCTUPOLE)
        .value("MAGNETIC_DIPOLE", OperatorType::MAGNETIC_DIPOLE)
        .value("DIAMAGNETIC", OperatorType::DIAMAGNETIC)
        .value("ARBITRARY", OperatorType::ARBITRARY);
}
