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

#ifndef NUMPYUTILS_H
#define NUMPYUTILS_H

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <Python.h>
#include <numpy/ndarrayobject.h>

#include <iterator>
#include <stdexcept>
#include <cstring>
#include <numeric>
#include <tuple>

#include <Eigen/Core>
#include <Eigen/SparseCore>

#include "Traits.h"

namespace numpy {

    /** \brief Selector for view or copy
     *
     * Deduction for view and copy should be a template parameter for
     * maximum flexibility.  That is why we make a new type for that.
     */
    enum view_or_copy { view, copy };

    /** \brief Numpy array type
     *
     * Numpy arrays are just usual PyObjects.  However, it doesn't
     * look as nice in return statements as numpy::array.
     */
    typedef PyObject* array;

    namespace internal {

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

        /** \brief Specialization of py_type for float */
        template < > struct py_type < float >
        {
            /** \brief Numpy type identifier */
            static constexpr int type = NPY_FLOAT;
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

        /** \brief Perform sanity checks for convert
         *
         * To ensure that dimensions and types are consistent we
         * introduce this check function.
         *
         * \param[in] len      Length of the memory
         * \param[in] nd       Number of dimensions
         * \param[in] dims     Initializer list with lengths of dimensions
         *
         * \return The total length of the array
         */
        inline void check_array_sanity(int len, int nd, std::initializer_list<long> dims)
        {
            if ( len < 1 )
              throw std::out_of_range(
                "Trying to create a numpy array with zero or negative element count!");

            if ( nd != static_cast<int>(dims.size()) )
              throw std::out_of_range(
                "Dimension mismatch!");

            if ( len > std::accumulate( dims.begin(), dims.end(), int{1}, [](int a, int b) { return a*b; } ) )
              throw std::out_of_range(
                "Requested dimension is larger than data!");
        }

        /** \brief Create a Numpy view of an iterator (implementation)
         *
         * This creates a view of the data, i.e. the data still
         * belongs to the C++ end.  When the data is deallocated on
         * the C++ end trying to access it in Python will probably
         * cause a segfault.  In any case it is undefined behavior.
         *
         * \param[in] nd       Number of dimensions
         * \param[in] dim      Dimensions of the array
         * \param[in] dtype    Value type of array elements
         * \param[in] data     Pointer to array data
         *
         * \return PyObject* containing a Numpy array
         */
        template < view_or_copy v, bool const_tag, typename value_type >
        typename std::enable_if < v == numpy::view, PyObject * >::type
        convert_impl(int nd, npy_intp * dim, int dtype, void * data, int)
        {
            PyObject * ndarray = PyArray_New(&PyArray_Type, nd, dim, dtype, nullptr,
                                             data, 0, NPY_ARRAY_FARRAY, nullptr);

            if ( const_tag )
                PyArray_CLEARFLAGS(reinterpret_cast<PyArrayObject*>(ndarray),
                                   NPY_ARRAY_WRITEABLE);
            return ndarray;
        }

        /** \brief Create a Numpy copy of an iterator (implementation)
         *
         * This creates a copy of the data, i.e. the data belongs to
         * the Python end.  The array is marked with NPY_OWNDATA which
         * should tell the garbage collector to release this memory.
         *
         * \param[in] nd       Number of dimensions
         * \param[in] dim      Dimensions of the array
         * \param[in] dtype    Value type of array elements
         * \param[in] data     Pointer to array data
         * \param[in] len      Length of data
         *
         * \return PyObject* containing a Numpy array
         */
        template < view_or_copy v, bool const_tag, typename value_type >
        typename std::enable_if < v == numpy::copy, PyObject * >::type
        convert_impl(int nd, npy_intp * dim, int dtype, void * data, int len)
        {
            PyObject * ndarray = PyArray_New(&PyArray_Type, nd, dim, dtype, nullptr,
                                             nullptr, 0, NPY_ARRAY_FARRAY, nullptr);

            std::memcpy(PyArray_DATA(reinterpret_cast<PyArrayObject*>(ndarray)),
                        data, len*sizeof(value_type));
            PyArray_ENABLEFLAGS(reinterpret_cast<PyArrayObject*>(ndarray),
                                NPY_ARRAY_OWNDATA);
            return ndarray;
        }

        /** \brief Create a Numpy array from an iterator (implementation)
         *
         * This is the implementation for the conversion of an
         * arbitrary pointer between \p begin and \p end .
         *
         * \param[in] begin    Random access iterator pointing to the start
         * \param[in] end      Random access iterator pointing to the end
         * \param[in] nd       Number of dimensions
         * \param[in] dims     Initializer list with lengths of dimensions
         *
         * \return PyObject* containing a Numpy array
         */
        template < view_or_copy v, typename RAIter >
        PyObject * convert(RAIter begin, RAIter end, int nd, std::initializer_list<long> dims,
                             std::random_access_iterator_tag)
        {
            using value_type = typename std::iterator_traits<RAIter>::value_type;
            constexpr auto const_tag = traits::is_pointer_to_const < decltype(&(*begin)) >::value;

            int const len = std::distance(begin, end);

            check_array_sanity(len, nd, dims);

            npy_intp * dim = const_cast <
                typename traits::pointer_remove_const<decltype(&(*dims.begin()))>::type
                > ( &(*dims.begin()) );
            auto dtype = py_type<value_type>::type;
            void * data = const_cast <
                typename traits::pointer_remove_const<decltype(&(*begin))>::type
                > ( &(*begin) );

            return convert_impl<v, const_tag, value_type>(nd, dim, dtype, data, len);
        }
    }

    // View and copy of bare pointes

    /** \brief Create a Numpy array from an iterator (frontend)
     *
     * This is a proxy for internal::view_impl to make sure that the
     * iterator passed is random access.
     *
     * \param[in] begin    Iterator pointing to the start
     * \param[in] end      Iterator pointing to the end
     * \param[in] nd       Number of dimensions
     * \param[in] dims     Initializer list with lengths of dimensions
     *
     * \return PyObject* containing a Numpy array
     */
    template < view_or_copy v, typename Iter >
    PyObject * convert(Iter begin, Iter end, int nd, std::initializer_list<long> dims)
    {
        return internal::convert<v>(
          begin, end, nd, dims,
          typename std::iterator_traits<Iter>::iterator_category());
    }

    // Overloads for often used types

    /** \brief Create a Numpy array from a raw pointer
     *
     * Specialization of numpy::convert for a raw pointer.
     *
     * \param[in] begin    pointer to beginning
     * \param[in] end      pointer to end
     *
     * \return PyObject* containing a Numpy array
     */
    template < view_or_copy v, typename T >
    PyObject * convert(T * begin, T * end)
    {
        return numpy::convert<v>(begin, end, 1, {std::distance(begin,end)});
    }
    template < view_or_copy v, typename T >
    PyObject * convert(T const * begin, T const * end)
    {
        return numpy::convert<v>(begin, end, 1, {std::distance(begin,end)});
    }

    /** \brief Create a Numpy array from a Eigen::Matrix
     *
     * Specialization of numpy::convert for Eigen::Matrix
     *
     * \param[in] m    Eigen::Matrix&
     *
     * \return PyObject* containing a Numpy array
     */
    template < view_or_copy v, typename T >
    PyObject * convert(Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic>& m)
    {
        return numpy::convert<v>(m.data(), m.data()+m.size(), 2, {m.rows(), m.cols()});
    }
    template < view_or_copy v, typename T >
    PyObject * convert(Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> const& m)
    {
        return numpy::convert<v>(m.data(), m.data()+m.size(), 2, {m.rows(), m.cols()});
    }

    /** \brief Create a Numpy array from an Eigen::Vector
     *
     * Specialization of numpy::convert for Eigen::Vector
     *
     * \param[in] vec    Eigen::Vector&
     *
     * \return PyObject* containing a Numpy array
     */
    template < view_or_copy v, typename T >
    PyObject * convert(Eigen::Matrix<T,Eigen::Dynamic,1>& vec)
    {
        return numpy::convert<v>(vec.data(), vec.data()+vec.size(), 1, {vec.size()});
    }
    template < view_or_copy v, typename T >
    PyObject * convert(Eigen::Matrix<T,Eigen::Dynamic,1> const& vec)
    {
        return numpy::convert<v>(vec.data(), vec.data()+vec.size(), 1, {vec.size()});
    }

    // Experimental support for Sparse matrix return types

    namespace internal {

        /** \brief Create a Numpy array from an Eigen::SparseMatrix
         *
         * This is the implementation for the conversion of an
         * Eigen::SparseMatrix to a scipy.sparse.csc_matrix.  We load
         * the scipy.sparse module and call the constructor for the
         * csc_matrix which we populate with numpy arrays.
         *
         * \param[in] sm    T&&
         *
         * \return PyObject* containing a Scipy csc_matrix
         */
        template < view_or_copy v, typename T >
        PyObject * sparse_impl(T&& sm)
        {
            if ( ! sm.isCompressed() )
                // we cannot call makeCompressed here because sm might be const
                throw std::runtime_error("Sparse matrix is not compressed!");

            numpy::array indptr = numpy::convert<v>(
                sm.outerIndexPtr(), sm.outerIndexPtr() + sm.outerSize() + 1);

            numpy::array indices = numpy::convert<v>(
                sm.innerIndexPtr(), sm.innerIndexPtr() + sm.nonZeros());

            numpy::array data = numpy::convert<v>(
                sm.valuePtr(), sm.valuePtr() + sm.nonZeros());

            int num_rows = sm.rows();
            int num_cols = sm.cols();

            char object[] = "csc_matrix";
            char arglist[] = "(OOO)(ii)";
            PyObject * scipy = PyImport_ImportModule("scipy.sparse");
            PyObject * mat = PyObject_CallMethod(
                scipy, object, arglist, data, indices, indptr, num_rows, num_cols);
            Py_DECREF(scipy);
            return mat;
        }

    } // namespace internal

    /** \brief Create a Numpy array from an Eigen::SparseMatrix
     *
     * Specialization of numpy::convert for Eigen::SparseMatrix
     *
     * \param[in] sm    Eigen::SparseMatrix&
     *
     * \return PyObject* containing a Numpy array
     */
    template < view_or_copy v, typename T >
    PyObject * convert(Eigen::SparseMatrix<T>& sm)
    {
        return internal::sparse_impl<v>(sm);
    }
    template < view_or_copy v, typename T >
    PyObject * convert(Eigen::SparseMatrix<T> const& sm)
    {
        return internal::sparse_impl<v>(sm);
    }

    // Experimental support for Numpy arguments

    namespace internal {

        /** \brief Check if a Numpy has given storage order and alignment
         *
         * To be convertible to an Eigen matrix type the memory of a
         * Numpy array has to be FORTRAN style contiguous and has to
         * be aligned.
         *
         * \param[in] ndarray    Numpy array
         */
        void check_order_and_alignment(PyArrayObject * ndarray)
        {
            int const flags = PyArray_FLAGS(ndarray);
            if ( (flags & NPY_ARRAY_F_CONTIGUOUS) == 0 )
                throw std::invalid_argument(
                    "The argument is not contiguous or has wrong storage order!");
            if ( (flags & NPY_ARRAY_ALIGNED) == 0 )
                throw std::invalid_argument(
                    "The argument is not not aligned!");
        }

        /** \brief Deconstruct a numpy array for conversion to dense Eigen type
         *
         * \param[in] ndarray_    Numpy array
         *
         * \return A tuple of the number of dimensions, an array of
         * dimensions, and a pointer to the data.
         */
        template < typename T, typename Scalar = typename T::Scalar >
        std::tuple<int, npy_intp *, Scalar *> as_dense_impl(numpy::array ndarray_)
        {
            if ( ! ndarray_ || ! PyArray_Check(ndarray_) )
                throw std::invalid_argument(
                    "The argument is not a valid Numpy array!");

            PyArrayObject* ndarray = reinterpret_cast<PyArrayObject*>(ndarray_);
            internal::check_order_and_alignment(ndarray);
            int const nd = PyArray_NDIM(ndarray);
            npy_intp * dims = PyArray_SHAPE(ndarray);
            Scalar * data = reinterpret_cast<Scalar *>(PyArray_DATA(ndarray));

            if ( PyArray_TYPE(ndarray) != internal::py_type<Scalar>::type )
                throw std::invalid_argument("Type mismatch.");

            return std::make_tuple(nd, dims, data);
        }

    } // namespace internal

    /** \brief Convert a Numpy array to an Eigen::Matrix
     *
     * \warning Experimental!  For now the returned Eigen::Map is
     * read-only.
     *
     * \param[in] ndarray    Numpy array
     *
     * \return An Eigen::Map of the Numpy array memory
     */
#ifndef SCANNED_BY_DOXYGEN
    template < typename T, typename Scalar = typename T::Scalar >
    typename std::enable_if<
        std::is_same<
            T,
            Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic, T::Options,
                          T::MaxRowsAtCompileTime, T::MaxColsAtCompileTime>
            >::value,
        Eigen::Map < T const >
        >::type
#else
    template < typename MatrixType >
    Eigen::Map < MatrixType >
#endif
    as(numpy::array ndarray)
    {
        int nd; npy_intp * dims; Scalar * data;
        std::tie(nd, dims, data) = internal::as_dense_impl<T>(ndarray);

        if ( nd > 2 )
            throw std::range_error("Dimension mismatch.");

        return Eigen::Map < T const > (data, dims[0], ( nd == 1 ? 1 : dims[1] ));
    }


    /** \brief Convert a Numpy array to an Eigen::Vector
     *
     * \warning Experimental!  For now the returned Eigen::Map is
     * read-only.
     *
     * \param[in] ndarray    Numpy array
     *
     * \return An Eigen::Map of the Numpy array memory
     */
#ifndef SCANNED_BY_DOXYGEN
    template < typename T, typename Scalar = typename T::Scalar >
    typename std::enable_if<
        std::is_same<
            T,
            Eigen::Matrix<Scalar, Eigen::Dynamic, 1, T::Options,
                          T::MaxRowsAtCompileTime, T::MaxColsAtCompileTime>
            >::value,
        Eigen::Map < T const >
        >::type
#else
    template < typename VectorType >
    Eigen::Map < VectorType >
#endif
    as(numpy::array ndarray)
    {
        int nd; npy_intp * dims; Scalar * data;
        std::tie(nd, dims, data) = internal::as_dense_impl<T>(ndarray);

        if ( nd > 1 )
            throw std::range_error("Dimension mismatch.");

        return Eigen::Map < T const > (data, dims[0]);
    }

    /** \brief Convert a Numpy array to an Eigen::SparseMatrix
     *
     * \warning Experimental!  For now the returned Eigen::Map is
     * read-only.
     *
     * \param[in] ndarray    Numpy array
     *
     * \return An Eigen::Map of the Numpy array memory
     */
#ifndef SCANNED_BY_DOXYGEN
    template < typename T, typename Scalar = typename T::Scalar,
               typename StorageIndex = typename T::StorageIndex >
    typename std::enable_if<
        std::is_same<
            T,
            Eigen::SparseMatrix<Scalar, T::Options, StorageIndex>
            >::value,
        Eigen::Map < T const >
        >::type
#else
    template < typename SparseMatrixType >
    Eigen::Map < SparseMatrixType >
#endif
    as(numpy::array ndarray)
    {
        PyObject * scipy = PyImport_ImportModule("scipy.sparse");
        PyObject * csc_matrix = PyObject_GetAttrString(scipy, "csc_matrix");
        if ( ! ndarray || ! PyObject_IsInstance(ndarray, csc_matrix) )
            throw std::invalid_argument(
                "The argument is not a valid csc_matrix!");
        Py_DECREF(csc_matrix);
        Py_DECREF(scipy);

        PyObject * tmp;

        tmp = PyObject_GetAttrString(ndarray, "data");
        PyArrayObject * data_ = reinterpret_cast<PyArrayObject*>(tmp);
        Py_DECREF(tmp); // Probably unsafe but otherwise the passed object can never be garbage collected
        internal::check_order_and_alignment(data_);

        tmp = PyObject_GetAttrString(ndarray, "indices");
        PyArrayObject * indices_ = reinterpret_cast<PyArrayObject*>(tmp);
        Py_DECREF(tmp); // Probably unsafe but otherwise the passed object can never be garbage collected
        internal::check_order_and_alignment(indices_);

        tmp = PyObject_GetAttrString(ndarray, "indptr");
        PyArrayObject * indptr_  = reinterpret_cast<PyArrayObject*>(tmp);
        Py_DECREF(tmp); // Probably unsafe but otherwise the passed object can never be garbage collected
        internal::check_order_and_alignment(indptr_);

        PyObject * shape = PyObject_GetAttrString(ndarray, "shape");
        int rows, cols;
        PyArg_ParseTuple(shape, "ii", &rows, &cols);
        Py_DECREF(shape);

        npy_intp * nnz = PyArray_SHAPE(data_);

        Scalar       * data    = reinterpret_cast<Scalar       *>(PyArray_DATA(data_   ));
        StorageIndex * indices = reinterpret_cast<StorageIndex *>(PyArray_DATA(indices_));
        StorageIndex * indptr  = reinterpret_cast<StorageIndex *>(PyArray_DATA(indptr_ ));

        if ( PyArray_TYPE(data_) != internal::py_type<Scalar>::type )
            throw std::invalid_argument("Type mismatch.");

        return Eigen::Map < T const > (rows, cols, *nnz, indptr, indices, data);
    }

}

#endif // NUMPYUTILS_H
