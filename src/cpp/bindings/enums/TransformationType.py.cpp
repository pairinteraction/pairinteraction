#include "TransformationType.py.hpp"

#include "TransformationType.hpp"

#include <nanobind/nanobind.h>

namespace nb = nanobind;

void bind_transformation_type(nb::module_ &m) {
    nb::enum_<TransformationType>(m, "TransformationType", nb::is_arithmetic())
        .value("NONE", TransformationType::NONE)
        .value("SORT_BY_KET", TransformationType::SORT_BY_KET)
        .value("SORT_BY_QUANTUM_NUMBER_F", TransformationType::SORT_BY_QUANTUM_NUMBER_F)
        .value("SORT_BY_QUANTUM_NUMBER_M", TransformationType::SORT_BY_QUANTUM_NUMBER_M)
        .value("SORT_BY_PARITY", TransformationType::SORT_BY_PARITY)
        .value("SORT_BY_ENERGY", TransformationType::SORT_BY_ENERGY)
        .value("ROTATE", TransformationType::ROTATE)
        .value("ARBITRARY", TransformationType::ARBITRARY)
        .value("MASK_SORTING", TransformationType::MASK_SORTING);
}
