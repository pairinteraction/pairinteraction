// SPDX-FileCopyrightText: 2024 PairInteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include <cstddef>

namespace pairinteraction::utils {
constexpr std::size_t M_OFFSET = 500;
constexpr std::size_t KET_ID_STRIDE = 2 * M_OFFSET;

struct IdAndM {
    std::size_t id;
    double m;
};

inline std::size_t encode_as_ket_id(IdAndM id_and_m) {
    return id_and_m.id * KET_ID_STRIDE + static_cast<std::size_t>(2 * id_and_m.m + M_OFFSET);
}

inline IdAndM decode_from_ket_id(std::size_t ket_id) {
    return {
        .id = ket_id / KET_ID_STRIDE,
        .m = (static_cast<double>(ket_id % KET_ID_STRIDE) - static_cast<double>(M_OFFSET)) / 2.0,
    };
}
} // namespace pairinteraction::utils
