/*
 * Copyright (c) 2017 Sebastian Weber, Henri Menke. All rights reserved.
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

#ifndef TRAITS_H
#define TRAITS_H

namespace traits {

    /** \brief Check if something points to const
     *
     * This struct has a member variable indicating whether a type
     * points to const
     */
    template < typename T >
    struct is_pointer_to_const : std::false_type {};

    /** \brief Specialization of is_pointer_to_const for T const * */
    template < typename T >
    struct is_pointer_to_const < T const * > : std::true_type {};

    /** \brief Specialization of is_pointer_to_const for T const * const */
    template < typename T >
    struct is_pointer_to_const < T const * const > : std::true_type {};

    /** \brief Remove const qualifiers from pointer
     *
     * This struct has a member typedef which holds the pointer
     * type without const qualifiers.
     */
    template < typename T >
    struct pointer_remove_const;

    /** \brief Specialization of pointer_remove_const for T* */
    template < typename T >
    struct pointer_remove_const < T * >
    {
        /** \brief Pointer type without const */
        typedef T* type;
    };

    /** \brief Specialization of pointer_remove_const for T const * */
    template<typename T>
    struct pointer_remove_const < T const * >
    {
        /** \brief Pointer type without const */
        typedef T* type;
    };

    /** \brief Specialization of pointer_remove_const for T const * const */
    template<typename T>
    struct pointer_remove_const < T const * const >
    {
        /** \brief Pointer type without const */
        typedef T* type;
    };

    /** \brief Add const qualifiers to pointer
     *
     * This struct has a member typedef which holds the pointer
     * type with additional const qualifiers.
     */
    template < typename T >
    struct pointer_add_const;

    /** \brief Specialization of pointer_add_const for T* */
    template < typename T >
    struct pointer_add_const < T * >
    {
        /** \brief Pointer type without const */
        typedef T const * const type;
    };

    /** \brief Specialization of pointer_add_const for T const * */
    template<typename T>
    struct pointer_add_const < T const * >
    {
        /** \brief Pointer type without const */
        typedef T const * const type;
    };

    /** \brief Specialization of pointer_add_const for T const * const */
    template<typename T>
    struct pointer_add_const < T const * const >
    {
        /** \brief Pointer type without const */
        typedef T const * const type;
    };

} // namespace traits

#endif // TRAITS_H
