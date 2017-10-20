/*
 * Copyright (c) 2016 Sebastian Weber, Henri Menke. All rights reserved.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#ifndef ITER_H
#define ITER_H

template<class ContainerType, class ValueType>
class ConstIter {
public:
    ConstIter (const ContainerType* p_vec, int pos)
      : _pos( pos ) , _p_vec( p_vec ) { }

    bool operator!= (ConstIter const& other) const
    {
        return _pos != other._pos;
    }

    ValueType operator* () const
    {
        return _p_vec->get( _pos );
    }

    ConstIter& operator++ ()
    {
        ++_pos;
        return *this;
    }

private:
    int _pos;
    const ContainerType* _p_vec;
};

#endif
