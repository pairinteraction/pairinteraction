#include "./run_unit_tests.py.hpp"

#include "pairinteraction/tools/run_unit_tests.hpp"

#include <nanobind/nanobind.h>
#include <nanobind/stl/filesystem.h>

namespace nb = nanobind;
using namespace pairinteraction;

void bind_run_unit_tests(nb::module_ &m) {
    m.def("run_unit_tests",
          [](bool download_missing, bool wigner_in_memory, std::filesystem::path database_dir) {
              run_unit_tests(0, nullptr, download_missing, wigner_in_memory,
                             std::move(database_dir));
          });
}
