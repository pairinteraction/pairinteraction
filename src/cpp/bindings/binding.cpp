#include "Ket.py.hpp"
#include "pairinteraction/tools/setup.hpp"

#include <nanobind/nanobind.h>

NB_MODULE(backend, m) { // NOLINT
    pairinteraction::setup();
    bind_ket(m);
}
