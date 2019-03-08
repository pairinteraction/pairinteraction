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

#include "SystemBase.h"

// Explicit template specializations

template <>
template <>
StateOne SystemBase<StateOne>::createStateFromLabel(const std::string &label) const {
    return StateOne(label);
}

template <>
template <>
StateTwo SystemBase<StateTwo>::createStateFromLabel(const std::string &label) const {
    return StateTwo(std::array<std::string, 2>({{"0_" + label, "1_" + label}}));
}
