// -*- C -*-  (not really, but good for syntax highlighting)
%module(docstring="Python interface for pairinteraction") binding

%{
#define SWIG_FILE_WITH_INIT

#include "Constants.hpp"
#include "Interface.hpp"
#include "MatrixElementCache.hpp"
#include "PerturbativeInteraction.hpp"
#include "QuantumDefect.hpp"
#include "State.hpp"
#include "Symmetry.hpp"
#include "SystemOne.hpp"
#include "SystemTwo.hpp"
#include "Wavefunction.hpp"

#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>

#include <streambuf>
#include <sstream>
#include <string.h>

struct array_source : std::streambuf {
    array_source(char *buffer, size_t len) {
        this->setg(buffer, buffer, buffer + len);
    }
};
%}


%feature("autodoc", "");


%include "attribute.i"
%include "exception.i"
%include "std_string.i"
%include "std_vector.i"

#ifdef SWIGPYTHON
%{
#include "NumpyUtils.hpp"
#include "Traits.hpp"
%}

%init %{
import_array();
%}

%include "std_array.i"
%include "std_complex.i"
%include "std_set.i"
#endif

// Convert C++ exceptions to Python exceptions
// http://www.swig.org/Doc1.3/Library.html#Library_stl_exceptions
%exception {
  try {
    $action
  } catch (const std::exception& e) {
    SWIG_exception(SWIG_RuntimeError, e.what());
  }
}

// Make pickle work
// https://stackoverflow.com/questions/9310053/how-to-make-my-swig-extension-module-work-with-pickle
// remark: passing as std::stringstream does not work because swig calls the implicitly-deleted copy constructor of std::stringstream instead of the move constructor
#ifdef SWIGPYTHON
%define %boost_picklable(cls...)
  %extend cls {
    PyObject* __getstate__() {
      std::stringstream ss;
      boost::archive::binary_oarchive ar(ss);
      ar << *($self);
      return PyBytes_FromStringAndSize(ss.str().data(), ss.str().length());
    }

    void __setstate_internal(PyObject* const sState) {
      char *buffer;
      Py_ssize_t len;
      PyBytes_AsStringAndSize(sState, &buffer, &len);
      array_source asource(buffer, static_cast<size_t>(len));
      std::istream ss(&asource);
      boost::archive::binary_iarchive ar(ss);
      ar >> *($self);
    }
  }
%enddef
#else
%define %boost_picklable(cls...)
%enddef
#endif
// Wrap the << operator
// http://stackoverflow.com/questions/2548779/how-to-stringfy-a-swig-matrix-object-in-python
%define __STR__()
std::string __str__() {
  std::ostringstream out;
  out << *$self;
  return out.str();
}
%enddef

 // Instantiate often used STL templates
namespace std {
  %template(VectorInt) vector<int>;
  %template(VectorDouble) vector<double>;
  %template(VectorFloat) vector<float>;
  %template(VectorStateOne) vector<StateOne>;
  %template(VectorStateTwo) vector<StateTwo>;
  %template(VectorSizeT) vector<size_t>;
  %template(VectorComplexDouble) vector<std::complex<double>>;
};

// Make numpy wrappers
#ifdef SWIGPYTHON
%define MAKE_NUMPY_TYPEMAP_OUT(TYPE...)
  %typemap(out) TYPE        { $result = numpy::convert<numpy::copy>( $1); }
  %typemap(out) TYPE &      { $result = numpy::convert<numpy::view>(*$1); }
  %typemap(out) TYPE const&
  {
    $result = numpy::convert<numpy::view>(*const_cast<
      traits::pointer_add_const< decltype($1) >::type >($1));
  }
%enddef

%define MAKE_NUMPY_TYPEMAP_IN(TYPE...)
  %typemap(in) TYPE // SWIG omits the exception handler here
  {
    try {
      $1 = numpy::as<TYPE>($input);
    } catch(std::exception &e) {
      SWIG_exception(SWIG_RuntimeError, e.what());
    }
  }
%enddef

%define MAKE_NUMPY_TYPEMAP(TYPE...)
  MAKE_NUMPY_TYPEMAP_OUT(TYPE)
  MAKE_NUMPY_TYPEMAP_IN(TYPE)
%enddef

MAKE_NUMPY_TYPEMAP(Eigen::MatrixX<double>)
MAKE_NUMPY_TYPEMAP(Eigen::VectorX<double>)
MAKE_NUMPY_TYPEMAP(Eigen::SparseMatrix<double>)
MAKE_NUMPY_TYPEMAP(Eigen::SparseMatrix<std::complex<double>>)

 // Instantiate often used STL templates
namespace std {
  %template(ArrayBoolTwo) array<bool,2>;
  %template(ArrayStringTwo) array<string,2>;
  %template(ArrayIntTwo) array<int,2>;
  %template(ArrayFloatTwo) array<float,2>;
  %template(ArrayDoubleTwo) array<double,2>;
  %template(ArrayDoubleThree) array<double,3>;
  %template(ArrayVectorSizeTTwo) array<std::vector<size_t>,2>;
  %template(VectorArraySizeTTwo) vector<std::array<size_t, 2>>;
  %template(ArrayEigenVectorDoubleTwo) array<Eigen::VectorX<double>,2>;
  %template(SetInt) set<int>;
  %template(SetFloat) set<float>;
  %template(SetStateOne) set<StateOne>;
  %template(SetStateTwo) set<StateTwo>;
};

MAKE_NUMPY_TYPEMAP_OUT(std::vector<int>)
MAKE_NUMPY_TYPEMAP_OUT(std::vector<float>)
MAKE_NUMPY_TYPEMAP_OUT(std::vector<double>)
MAKE_NUMPY_TYPEMAP_OUT(std::array<int,2>)
MAKE_NUMPY_TYPEMAP_OUT(std::array<float,2>)
MAKE_NUMPY_TYPEMAP_OUT(std::array<double,2>)
#endif


%feature("autodoc", "2");


%rename(__lt__) operator<;

%include "Constants.hpp"
%include "Interface.hpp"
%include "Symmetry.hpp"

%template(computeComplex) compute<std::complex<double>>;
%template(computeReal) compute<double>;

%include "QuantumDefect.hpp"

%include "Wavefunction.hpp"

%include "PerturbativeInteraction.hpp"

%include "WignerD.hpp"

// Wrap MatrixElementCache.h
%include "MatrixElementCache.hpp"

%boost_picklable(MatrixElementCache);

%extend MatrixElementCache {
#ifdef SWIGPYTHON
  %pythoncode %{
    def __setstate__(self, sState):
      self.__init__()
      self.__setstate_internal(sState)
  %}
#endif
}


// Wrap State.h
%ignore hash;
%rename(__ostream__) operator<<;

%include "State.hpp"

%boost_picklable(StateOne);
%boost_picklable(StateTwo);

%extend StateOne {
  __STR__();
#ifdef SWIGPYTHON
  %pythoncode %{
    def __setstate__(self, sState):
      self.__init__()
      self.__setstate_internal(sState)
  %}
#endif
}


%extend StateTwo {
  __STR__();
#ifdef SWIGPYTHON
  %pythoncode %{
    def __setstate__(self, sState):
      self.__init__()
      self.__setstate_internal(sState)
  %}
#endif
}


// Wrap SystemBase.h
%warnfilter(509) SystemBase::getOverlap;

%include "SystemBase.hpp"


// Wrap SystemOne.h and SystemTwo.h
%template(_SystemStateOneComplex) SystemBase<std::complex<double>,StateOne>;
%template(_SystemStateTwoComplex) SystemBase<std::complex<double>,StateTwo>;
%template(_SystemStateOneReal) SystemBase<double,StateOne>;
%template(_SystemStateTwoReal) SystemBase<double,StateTwo>;

%copyctor SystemOne;
%copyctor SystemTwo;

%include "SystemOne.hpp"
%include "SystemTwo.hpp"

%boost_picklable(SystemOne);
%boost_picklable(SystemTwo);

%extend SystemOne {
#ifdef SWIGPYTHON
  %pythoncode %{
    def __setstate__(self, sState):
      tmp = MatrixElementCache()
      self.__init__("", tmp)
      self.__setstate_internal(sState)
  %}
#endif
}

%extend SystemTwo<double> {
#ifdef SWIGPYTHON
  %pythoncode %{
    def __setstate__(self, sState):
      tmp = MatrixElementCache()
      s1 = SystemOneReal("", tmp)
      s2 = SystemOneReal("", tmp)
      self.__init__(s1, s2, tmp)
      self.__setstate_internal(sState)
  %}
#endif
}

%extend SystemTwo<std::complex<double>> {
#ifdef SWIGPYTHON
  %pythoncode %{
    def __setstate__(self, sState):
      tmp = MatrixElementCache()
      s1 = SystemOneComplex("", tmp)
      s2 = SystemOneComplex("", tmp)
      self.__init__(s1, s2, tmp)
      self.__setstate_internal(sState)
  %}
#endif
}

%template(SystemOneComplex) SystemOne<std::complex<double>>;
%template(SystemTwoComplex) SystemTwo<std::complex<double>>;
%template(SystemOneReal) SystemOne<double>;
%template(SystemTwoReal) SystemTwo<double>;
