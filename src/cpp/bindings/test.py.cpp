#include "test.py.hpp"

#include "test.hpp"

#include <nanobind/nanobind.h>
#include <nanobind/stl/filesystem.h>

namespace nb = nanobind;

void bind_test(nb::module_ &m) {
    m.def("test", [](bool download_missing, std::filesystem::path databasedir) {
        test(0, nullptr, download_missing, databasedir);
    });
}
