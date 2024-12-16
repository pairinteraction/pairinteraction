#include "./version.py.hpp"

#include "pairinteraction/version.hpp"

namespace nb = nanobind;
using namespace pairinteraction;

void bind_version(nb::module_ &m) {
    m.attr("VERSION_MAJOR") = pairinteraction::VERSION_MAJOR;
    m.attr("VERSION_MINOR") = pairinteraction::VERSION_MINOR;
    m.attr("VERSION_PATCH") = pairinteraction::VERSION_PATCH;
}
