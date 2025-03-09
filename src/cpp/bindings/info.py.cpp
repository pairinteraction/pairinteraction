#include "pairinteraction/info.hpp"

#include "./info.py.hpp"

#include <nanobind/nanobind.h>

namespace nb = nanobind;
using namespace pairinteraction;

void bind_info(nb::module_ &m) {
    nb::class_<Info>(m, "Info")
        .def_ro("has_eigen", &Info::has_eigen)
        .def_ro("has_lapacke_evd", &Info::has_lapacke_evd)
        .def_ro("has_lapacke_evr", &Info::has_lapacke_evr)
        .def_ro("has_feast", &Info::has_feast);
}
