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
    MatrixCOO(size_t nRows, size_t nCols, std::vector<triple> &data);
    MatrixCOO();
    void add(didx row, didx col, dreal val);
    void multiplyScalar(dreal &&scalar);
    void order();
    void sumup();
    void print();
    std::shared_ptr<MatrixCRS> toCRS();
    std::vector<triple> vectorize();
    void devectorize(std::vector<triple> &vector);
    size_t getNumRows();
    size_t getNumCols();
    std::vector<size_t> getDimensions();
    void setDimensions(std::vector<size_t> &&dimensions);

private:
    size_t nRows;
    size_t nCols;
    std::vector<triple> data; // row, col, val
    bool ordered;
    bool sumuped;
};

#endif
