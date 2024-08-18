#include "pintr/ket/Ket.py.hpp"
#include "pintr/tools/setup.hpp"

#include <nanobind/nanobind.h>

NB_MODULE(backend, m) { // NOLINT
    pintr::setup();
    bind_ket(m);
}
