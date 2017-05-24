%include "std_complex.i"

%{
#include "Eigen/Core"
#include "Eigen/SparseCore"
%}

#define Index int

namespace Eigen {
    template<typename Scalar, int RowsAtCompileTime, int ColsAtCompileTime>
    class Matrix
    {
    public:
        Matrix(Index rows, Index cols);
    };

    template<typename Scalar>
    class SparseMatrix
    {
    public:
        SparseMatrix(Index rows, Index cols);
    };
}


%define DECLARE_COMMON_METHODS()
    Index rows() const { return self->rows(); };
    Index cols() const { return self->cols(); };
    std::string __str__() const {
        std::ostringstream out;
        out << *self;
        return out.str();
    }
%enddef


%template(VectorXd) Eigen::Matrix<double,Eigen::Dynamic,1>;
%extend Eigen::Matrix<double,Eigen::Dynamic,1> {
    DECLARE_COMMON_METHODS();
    double _get(Index i) const { return (*self)(i); }
    void _set(Index i, double x) { (*self)(i) = x; }
}

%template(VectorXcd) Eigen::Matrix<std::complex<double>,Eigen::Dynamic,1>;
%extend Eigen::Matrix<std::complex<double>,Eigen::Dynamic,1> {
    DECLARE_COMMON_METHODS();
    std::complex<double> _get(Index i) const { return (*self)(i); }
    void _set(Index i, std::complex<double> const& x) { (*self)(i) = x; }
}

%template(MatrixXd) Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic>;
%extend Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> {
    DECLARE_COMMON_METHODS();
    double _get(Index i, Index j) const { return (*self)(i,j); }
    void _set(Index i, Index j, double x) { (*self)(i,j) = x; }
}

%template(MatrixXcd) Eigen::Matrix<std::complex<double>,Eigen::Dynamic,Eigen::Dynamic>;
%extend Eigen::Matrix<std::complex<double>,Eigen::Dynamic,Eigen::Dynamic> {
    DECLARE_COMMON_METHODS();
    std::complex<double> _get(Index i, Index j) const { return (*self)(i,j); }
    void _set(Index i, Index j, std::complex<double> const& x) { (*self)(i,j) = x; }
}

%template(SparseMatrixXd) Eigen::SparseMatrix<double>;
%extend Eigen::SparseMatrix<double> {
    DECLARE_COMMON_METHODS();
    double _get(Index i, Index j) const { return self->coeff(i,j); }
    void _set(Index i, Index j, double x) { self->coeffRef(i,j) = x; }
}

%template(SparseMatrixXcd) Eigen::SparseMatrix<std::complex<double>>;
%extend Eigen::SparseMatrix<std::comlex<double>> {
    DECLARE_COMMON_METHODS();
    std::complex<double> _get(Index i, Index j) const { return self->coeff(i,j); }
    void _set(Index i, Index j, std::complex<double> const& x) { self->coeffRef(i,j) = x; }
}

// Python specific extensions

// This block has to be the last in the Interface.i file, because it
// will be copied literally into the Python module and we want to
// overwrite specific methods, defined before.

%pythoncode %{

def getitem_from_tuple(this,a):
    if a[0] >= this.rows() or a[1] >= this.cols():
        raise IndexError("index out of range")
    return this._get(a[0],a[1])

def setitem_from_tuple(this,a,x):
    if a[0] >= this.rows() or a[1] >= this.cols():
        raise IndexError("index out of range")
    this._set(a[0],a[1],x)

VectorXd .__getitem__ = VectorXd ._get
VectorXd .__setitem__ = VectorXd ._set
VectorXcd.__getitem__ = VectorXcd._get
VectorXcd.__setitem__ = VectorXcd._set
MatrixXd .__getitem__ = getitem_from_tuple
MatrixXd .__setitem__ = setitem_from_tuple
MatrixXcd.__getitem__ = getitem_from_tuple
MatrixXcd.__setitem__ = setitem_from_tuple
SparseMatrixXd .__getitem__ = getitem_from_tuple
SparseMatrixXd .__setitem__ = setitem_from_tuple
SparseMatrixXcd.__getitem__ = getitem_from_tuple
SparseMatrixXcd.__setitem__ = setitem_from_tuple

%}
