/*
 * Copyright (c) 2016 Sebastian Weber, Henri Menke. All rights reserved.
 *
 * This file is part of the pairinteraction library.
 *
 * The pairinteraction library is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * The pairinteraction library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with the pairinteraction library. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef SYMMETRY_H
#define SYMMETRY_H

#include <array>
#include <cstddef>
#include <climits>

enum parity_t {
    NA = INT_MAX,
    EVEN = 1,
    ODD = -1,
};

struct Symmetry {
    parity_t inversion;
    parity_t reflection;
    parity_t permutation;
    int rotation;

    // Comparison operator that is needed if an object of type Symmetry is used as key for std::map
    friend bool operator<(const Symmetry &s1, const Symmetry &s2) {
        std::array<int, 5> syms1{{s1.inversion, s1.reflection, s1.permutation, s1.rotation}};
        std::array<int, 5> syms2{{s2.inversion, s2.reflection, s2.permutation, s2.rotation}};

        for (size_t i = 0; i < syms1.size(); ++i) {
            if (syms1[i] < syms2[i]) {
                return true;
            }
            if (syms1[i] > syms2[i]) {
                return false;
            }
        }
        return false;
    }
};

#endif
