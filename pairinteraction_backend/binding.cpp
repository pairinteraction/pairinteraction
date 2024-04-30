#include "ket/Ket.py.hpp"

#include <nanobind/nanobind.h>

NB_MODULE(pairinteraction_backend, m) { bind_ket(m); }
