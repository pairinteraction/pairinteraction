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

#ifndef ITER_H
#define ITER_H

template <class ContainerType, class ValueType>
class ConstIter {
public:
    ConstIter(const ContainerType *p_vec, int pos) : _pos(pos), _p_vec(p_vec) {}

    bool operator!=(ConstIter const &other) const { return _pos != other._pos; }

    ValueType operator*() const { return _p_vec->get(_pos); }

    ConstIter &operator++() {
        ++_pos;
        return *this;
    }

private:
    int _pos;
    const ContainerType *_p_vec;
};

#endif
