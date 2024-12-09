#include "./test.py.hpp"

#include "pairinteraction/tools/test.hpp"

#include <nanobind/nanobind.h>
#include <nanobind/stl/filesystem.h>

namespace nb = nanobind;
using namespace pairinteraction;

void bind_test(nb::module_ &m) {
    m.def("test",
          [](bool download_missing, bool wigner_in_memory, std::filesystem::path database_dir) {
              test(0, nullptr, download_missing, wigner_in_memory, std::move(database_dir));
          });
}
