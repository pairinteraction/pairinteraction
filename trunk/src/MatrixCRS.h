#ifndef MATRIXCRS_H
#define MATRIXCRS_H

#include "dtypes.h"
#include "Serializable.h"
#include "Iter.h"

#include <math.h>
#include <vector>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <memory>
#include <unordered_map>

class MatrixCOO;

class RowCRS {
public:
    RowCRS (idx_t row, std::vector<idx_t> cols, std::vector<real_t> vals) {
        row_ = row;
        col_ = cols;
        val_ = vals;
    }
    Triple get (idx_t i) const {
        return Triple(row_,col_[i],val_[i]);
    }
    Iter<RowCRS, Triple> begin () const {
        return Iter<RowCRS, Triple>( this, 0 );
    }
    Iter<RowCRS, Triple> end () const {
        return Iter<RowCRS, Triple>( this, val_.size() );
    }
    void set (idx_t i, int v) {
        val_[i] = v;
    }

private:
    idx_t row_;
    std::vector<idx_t> col_;
    std::vector<real_t> val_;
};

class MatrixCRS : public Serializable {
public:
    MatrixCRS(size_t nRows, size_t nCols, size_t size);
    MatrixCRS(bytes_t &bytes);
    MatrixCRS();
    void add(idx_t r, idx_t c, real_t v);
    void multiplyScalar(real_t &&scalar);
    void order();
    void sumup();
    void print();
    void init(size_t nRows, size_t nCols, size_t size);
    RowCRS row(size_t r);
    MatrixCRS transformBackward(MatrixCRS &basis);
    MatrixCRS transformForward(MatrixCRS &basisoriginal);
    MatrixCOO toCOO();
    bytes_t serialize();
    void deserialize(bytes_t &bytes);
    size_t getNumRows();
    size_t getNumCols();
    size_t size();
    std::vector<size_t> getDimensions();
    void setDimensions(std::vector<size_t> &&dimensions);

    Triple get (idx_t i) const {
        return Triple(row_[i],col_[i],val_[i]);
    }
    Iter<MatrixCRS, Triple> begin () const {
        return Iter<MatrixCRS, Triple>( this, 0 );
    }
    Iter<MatrixCRS, Triple> end () const {
        return Iter<MatrixCRS, Triple>( this, val_.size() );
    }
    void set (idx_t i, int v) {
        val_[i] = v;
    }

private:
    void buildRow();

    std::vector<idx_t> ptr_;
    std::vector<idx_t> col_;
    std::vector<idx_t> row_;
    std::vector<real_t> val_;
    size_t nRows;
    size_t nCols;
    size_t r_last;
    bool ordered;
    bool sumuped;
};

#endif
