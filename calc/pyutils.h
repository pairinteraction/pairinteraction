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

#ifndef PYUTILS_H
#define PYUTILS_H

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <Python.h>
#include <numpy/ndarrayobject.h>
#include <iterator>
#include <stdexcept>
#include <cstring>
#include <vector>
#include <Eigen/Core>

namespace numpy {

    /** \brief Numpy array type
     *
     * Numpy arrays are just usual PyObjects.  However, it doesn't
     * look as nice in return statements as numpy::array.
     */
    typedef PyObject* array;

    namespace internal {

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

        /** \brief Map C++ types to Numpy types
         *
         * This struct has specializations for common C++ types and
         * their Numpy counterparts.
         */
        template < typename T >
        struct py_type;

        /** \brief Specialization of py_type for int */
        template < > struct py_type < int >
        {
            /** \brief Numpy type identifier */
            static constexpr int type = NPY_INT;
        };

        /** \brief Specialization of py_type for double */
        template < > struct py_type < double >
        {
            /** \brief Numpy type identifier */
            static constexpr int type = NPY_DOUBLE;
        };

        /** \brief Specialization of py_type for double */
        template < > struct py_type < std::complex<double> >
        {
            /** \brief Numpy type identifier */
            static constexpr int type = NPY_CDOUBLE;
        };

        /** \brief Create a Numpy view of an iterator (implementation)
         *
         * This is the implementation for the conversion of an
         * arbitrary pointer between \p begin and \p end with the
         * stride \p stride .  This creates a view of the data,
         * i.e. the data still belongs to the C++ end.  When the data
         * is deallocated on the C++ end trying to access it in Python
         * will probably cause a segfault.  In any case it is
         * undefined behavior.
         *
         * \param[in] begin    Random access iterator pointing to the start
         * \param[in] end      Random access iterator pointing to the end
         * \param[in] stride   Stride between rows of the data
         *
         * \return PyObject* containing a Numpy array
         */
        template < typename RAIter >
        PyObject * view_impl(RAIter begin, RAIter end, unsigned stride,
                             std::random_access_iterator_tag)
        {
            using value_type = typename std::iterator_traits<RAIter>::value_type;

            int const len  = std::distance(begin, end);
            if ( len < 1 )
              throw std::out_of_range(
                "Trying to create a numpy array with zero or negative element count");

            int const rows = len / stride;
            int const nd   = ( rows == 1 ? 1 : 2 );

            npy_intp dim[2] = { rows, stride };
            if ( nd == 1 ) dim[0] = { stride };
            auto dtype = internal::py_type<value_type>::type;
            void * data  = &(*begin);

            return PyArray_New(&PyArray_Type, nd, dim, dtype, nullptr,
                               data, 0, NPY_ARRAY_FARRAY, nullptr);
        }

        /** \brief Create a Numpy copy of an iterator (implementation)
         *
         * This is the implementation for the conversion of an
         * arbitrary pointer between \p begin and \p end with the
         * stride \p stride .  This creates a copy of the data,
         * i.e. the data belongs to the Python end.  The array is
         * marked with NPY_OWNDATA which should tell the garbage
         * collector to release this memory.
         *
         * \param[in] begin    Random access iterator pointing to the start
         * \param[in] end      Random access iterator pointing to the end
         * \param[in] stride   Stride between rows of the data
         *
         * \return PyObject* containing a Numpy array
         */
        template < typename RAIter >
        PyObject * copy_impl(RAIter begin, RAIter end, unsigned stride,
                             std::random_access_iterator_tag)
        {
            using value_type = typename std::iterator_traits<RAIter>::value_type;

            int const len  = std::distance(begin, end);
            if ( len < 1 )
              throw std::out_of_range(
                "Trying to create a numpy array with zero or negative element count");

            int const rows = len / stride;
            int const nd   = ( rows == 1 ? 1 : 2 );

            npy_intp dim[2] = { rows, stride };
            if ( nd == 1 ) dim[0] = { stride };
            auto dtype = internal::py_type<value_type>::type;
            void * data = const_cast <
                typename pointer_remove_const<decltype(&(*begin))>::type
                > ( &(*begin) );

            PyObject * ndarray = PyArray_New(&PyArray_Type, nd, dim, dtype, nullptr,
                                             nullptr, 0, NPY_ARRAY_FARRAY, nullptr);

            std::memcpy(PyArray_DATA(reinterpret_cast<PyArrayObject*>(ndarray)),
                        data, len*sizeof(value_type));
            PyArray_ENABLEFLAGS(reinterpret_cast<PyArrayObject*>(ndarray),
                                NPY_ARRAY_OWNDATA);
            return ndarray;
        }
    }

    // View and copy of bare pointes

    /** \brief Create a Numpy view of an iterator (frontend)
     *
     * This is a proxy for internal::view_impl to make sure that the
     * iterator passed is random access.
     *
     * \param[in] begin    Iterator pointing to the start
     * \param[in] end      Iterator pointing to the end
     * \param[in] stride   Stride between rows of the data
     *
     * \return PyObject* containing a Numpy array
     */
    template < typename Iter >
    PyObject * view(Iter begin, Iter end, unsigned stride)
    {
        return internal::view_impl(
          begin, end, stride,
          typename std::iterator_traits<Iter>::iterator_category());
    }

    /** \brief Create a Numpy copy of an iterator (frontend)
     *
     * This is a proxy for internal::copy_impl to make sure that the
     * iterator passed is random access.
     *
     * \param[in] begin    Iterator pointing to the start
     * \param[in] end      Iterator pointing to the end
     * \param[in] stride   Stride between rows of the data
     *
     * \return PyObject* containing a Numpy array
     */
    template < typename Iter >
    PyObject * copy(Iter begin, Iter end, unsigned stride)
    {
        return internal::copy_impl(
          begin, end, stride,
          typename std::iterator_traits<Iter>::iterator_category());
    }

    // Overloads for often used types

    /** \brief Create a Numpy view of a raw pointer
     *
     * Specialization of numpy::view for a raw pointer.
     *
     * \param[in] begin    pointer to beginning
     * \param[in] end      pointer to end
     *
     * \return PyObject* containing a Numpy array
     */
    template < typename T >
    PyObject * view(T * begin, T * end)
    {
        return numpy::view(begin, end, std::distance(begin,end));
    }

    /** \brief Create a Numpy copy of a raw pointer
     *
     * Specialization of numpy::view for a raw pointer.
     *
     * \param[in] begin    pointer to beginning
     * \param[in] end      pointer to end
     *
     * \return PyObject* containing a Numpy array
     */
    template < typename T >
    PyObject * copy(T * begin, T * end)
    {
        return numpy::copy(begin, end, std::distance(begin,end));
    }

    /** \brief Create a Numpy view of a std::vector
     *
     * Specialization of numpy::view for std::vector.
     *
     * \param[in] v    std::vector&
     *
     * \return PyObject* containing a Numpy array
     */
    template < typename T >
    PyObject * view(std::vector<T>& v)
    {
        return numpy::view(v.begin(), v.end(), v.size());
    }

    /** \brief Create a Numpy copy of a std::vector
     *
     * Specialization of numpy::view for std::vector.
     *
     * \param[in] v    std::vector&
     *
     * \return PyObject* containing a Numpy array
     */
    template < typename T >
    PyObject * copy(std::vector<T>& v)
    {
        return numpy::copy(v.begin(), v.end(), v.size());
    }

    /** \brief Create a Numpy view of a Eigen::Matrix
     *
     * Specialization of numpy::view for Eigen::Matrix
     *
     * \param[in] m    Eigen::Matrix&
     *
     * \return PyObject* containing a Numpy array
     */
    template < typename T >
    PyObject * view(Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic>& m)
    {
        return numpy::view(m.data(), m.data()+m.size(), m.cols());
    }

    /** \brief Create a Numpy copy of an Eigen::Matrix
     *
     * Specialization of numpy::view for Eigen::Matrix
     *
     * \param[in] m    Eigen::Matrix&
     *
     * \return PyObject* containing a Numpy array
     */
    template < typename T >
    PyObject * copy(Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic>& m)
    {
        return numpy::copy(m.data(), m.data()+m.size(), m.cols());
    }

    /** \brief Create a Numpy view of an Eigen::Vector
     *
     * Specialization of numpy::view for Eigen::Vector
     *
     * \param[in] v    Eigen::Vector&
     *
     * \return PyObject* containing a Numpy array
     */
    template < typename T >
    PyObject * view(Eigen::Matrix<T,Eigen::Dynamic,1>& v)
    {
        return numpy::view(v.data(), v.data()+v.size(), v.size());
    }

    /** \brief Create a Numpy copy of an Eigen::Vector
     *
     * Specialization of numpy::view for Eigen::Vector
     *
     * \param[in] v    Eigen::Vector&
     *
     * \return PyObject* containing a Numpy array
     */
    template < typename T >
    PyObject * copy(Eigen::Matrix<T,Eigen::Dynamic,1>& v)
    {
        return numpy::copy(v.data(), v.data()+v.size(), v.size());
    }
}

#endif // PYUTILS_H
