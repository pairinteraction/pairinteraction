#ifndef MATRIXCRS_H
#define MATRIXCRS_H

#include "dtypes.h"
#include "Serializable.h"

#include <math.h>
#include <vector>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <memory>

class MatrixCOO;

class MatrixCRS : public Serializable {
public:
    MatrixCRS(size_t nRows, size_t nCols, size_t size);
    MatrixCRS(bytes_t &bytes);
    MatrixCRS();
    void add(idx_t rIncrement, idx_t c, real_t v);
    void multiplyScalar(real_t &&scalar);
    void order();
    void sumup();
    void print();
    std::shared_ptr<MatrixCOO> toCOO();
    bytes_t serialize();
    void deserialize(bytes_t &bytes);
    size_t getNumRows();
    size_t getNumCols();
    std::vector<size_t> getDimensions();
    void setDimensions(std::vector<size_t> &&dimensions);

private:
    std::vector<idx_t> ptr;
    std::vector<idx_t> col;
    std::vector<real_t> val;
    size_t nRows;
    size_t nCols;
    bool ordered;
    bool sumuped;
};

#endif
