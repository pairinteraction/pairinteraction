#include "ket/Ket.py.hpp"
#include "setup.hpp"

#include <nanobind/nanobind.h>

NB_MODULE(pairinteraction_backend, m) { // NOLINT
    setup();
    bind_ket(m);
}
