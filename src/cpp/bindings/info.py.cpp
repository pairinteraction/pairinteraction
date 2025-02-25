#include "pairinteraction/info.hpp"

#include "./info.py.hpp"

namespace nb = nanobind;
using namespace pairinteraction;

void bind_info(nb::module_ &m) {
    nb::class_<Info>(m, "Info")
        .def_prop_ro_static("has_eigen", [](nb::object) { return Info::has_eigen; })
        .def_prop_ro_static("has_lapacke_evd", [](nb::object) { return Info::has_lapacke_evd; })
        .def_prop_ro_static("has_lapacke_evr", [](nb::object) { return Info::has_lapacke_evr; })
        .def_prop_ro_static("has_feast", [](nb::object) { return Info::has_feast; });
}
