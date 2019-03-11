/*
 * Copyright (c) 2017 Sebastian Weber, Henri Menke. All rights reserved.
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

#ifndef TRAITS_H
#define TRAITS_H

namespace traits {

/** \brief Check if something points to const
 *
 * This struct has a member variable indicating whether a type
 * points to const
 */
template <typename T>
struct is_pointer_to_const : std::false_type {};

/** \brief Specialization of is_pointer_to_const for T const * */
template <typename T>
struct is_pointer_to_const<T const *> : std::true_type {};

/** \brief Specialization of is_pointer_to_const for T const * const */
template <typename T>
struct is_pointer_to_const<T const *const> : std::true_type {};

/** \brief Remove const qualifiers from pointer
 *
 * This struct has a member typedef which holds the pointer
 * type without const qualifiers.
 */
template <typename T>
struct pointer_remove_const;

/** \brief Specialization of pointer_remove_const for T* */
template <typename T>
struct pointer_remove_const<T *> {
    /** \brief Pointer type without const */
    typedef T *type;
};

/** \brief Specialization of pointer_remove_const for T const * */
template <typename T>
struct pointer_remove_const<T const *> {
    /** \brief Pointer type without const */
    typedef T *type;
};

/** \brief Specialization of pointer_remove_const for T const * const */
template <typename T>
struct pointer_remove_const<T const *const> {
    /** \brief Pointer type without const */
    typedef T *type;
};

/** \brief Add const qualifiers to pointer
 *
 * This struct has a member typedef which holds the pointer
 * type with additional const qualifiers.
 */
template <typename T>
struct pointer_add_const;

/** \brief Specialization of pointer_add_const for T* */
template <typename T>
struct pointer_add_const<T *> {
    /** \brief Pointer type without const */
    typedef T const *const type;
};

/** \brief Specialization of pointer_add_const for T const * */
template <typename T>
struct pointer_add_const<T const *> {
    /** \brief Pointer type without const */
    typedef T const *const type;
};

/** \brief Specialization of pointer_add_const for T const * const */
template <typename T>
struct pointer_add_const<T const *const> {
    /** \brief Pointer type without const */
    typedef T const *const type;
};

} // namespace traits

#endif // TRAITS_H
