#ifndef MATRIXCOO_H
#define MATRIXCOO_H

#include "dtypes.h"
#include "Vectorizable.h"

#include <math.h>
#include <vector>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <memory>

class MatrixCRS;

class MatrixCOO : public Vectorizable {
public:
    MatrixCOO(size_t nRows, size_t nCols, size_t size);
    MatrixCOO(size_t nRows, size_t nCols, std::vector<Triple> &data);
    MatrixCOO();
    void add(idx_t row, idx_t col, real_t val);
    void multiplyScalar(real_t &&scalar);
    void order();
    void sumup();
    void print();
    std::shared_ptr<MatrixCRS> toCRS();
    std::vector<Triple> vectorize();
    void devectorize(std::vector<Triple> &vector);
    size_t getNumRows();
    size_t getNumCols();
    std::vector<size_t> getDimensions();
    void setDimensions(std::vector<size_t> &&dimensions);

private:
    size_t nRows;
    size_t nCols;
    std::vector<Triple> data; // row, col, val
    bool ordered;
    bool sumuped;
};

#endif
