#include "./Parity.py.hpp"

#include "pairinteraction/enums/Parity.hpp"

#include <nanobind/nanobind.h>

namespace nb = nanobind;
using namespace pairinteraction;

void bind_parity(nb::module_ &m) {
    nb::enum_<Parity>(m, "Parity")
        .value("ODD", Parity::ODD)
        .value("EVEN", Parity::EVEN)
        .value("UNKNOWN", Parity::UNKNOWN);
}
