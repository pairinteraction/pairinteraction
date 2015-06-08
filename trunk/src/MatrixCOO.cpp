#include "MatrixCOO.h"
#include "MatrixCRS.h"

MatrixCOO::MatrixCOO(size_t nRows, size_t nCols, size_t size) :  nRows(nRows), nCols(nCols), ordered(true), sumuped(true) { // empty
    data.reserve(size);
}

MatrixCOO::MatrixCOO(size_t nRows, size_t nCols, std::vector<Triple> &data) :  nRows(nRows), nCols(nCols) { // from vector
    devectorize(data);
}

MatrixCOO::MatrixCOO() {
}

void MatrixCOO::add(didx row, didx col, dreal val) {
    data.push_back(Triple(row,col,val));
    ordered = false;
    sumuped = false;
}


void MatrixCOO::multiplyScalar(dreal &&scalar) {
    for (auto &d: data) {
        d.val *= scalar;
    }
}

void MatrixCOO::order() {
    if (!ordered) {
        // Sort everything
        std::sort(data.begin(), data.end(),[](Triple t1, Triple t2) {return t1.col < t2.col;});
        std::stable_sort(data.begin(), data.end(),[](Triple t1, Triple t2) {return t1.row < t2.row;});

        ordered = true;
    }
}

void MatrixCOO::sumup() {
    order();

    if (!sumuped) {

        // Copy
        auto tripleTmp(data);
        data.clear();
        data.reserve(tripleTmp.size());

        // Sum doubled entries
        didx oldRow = 0;
        didx oldCol = 0;
        dreal sumVal = 0;
        for (auto &trp: tripleTmp) {
            didx row = trp.row;
            didx col = trp.col;
            dreal val = trp.val;

            if (row == oldRow && col == oldCol) {
                sumVal += val;
            } else {
                if (std::abs(sumVal) > 1e-12) {
                    data.push_back(Triple(oldRow,oldCol,sumVal));
                }
                sumVal = val;
                oldRow = row;
                oldCol = col;
            }
        }
        if (std::abs(sumVal) > 1e-12) {
            data.push_back(Triple(oldRow,oldCol,sumVal));
        }

        sumuped = true;
    }
}

void MatrixCOO::print() {
    order();
    sumup();

    didx n = 0;

    didx row = data[n].row;
    didx col = data[n].col;
    dreal val = data[n].val;

    std::cout << std::setiosflags(std::ios::right) << std::setiosflags(std::ios::fixed);
    for (didx r = 0; r < nRows; ++r) {
        for (didx c = 0; c < nCols; ++c) {
            dreal v = 0;
            if (row == r && col == c) {
                v = val;
                row = data[++n].row;
                col = data[n].col;
                val = data[n].val;
            }
            std::cout << std::setprecision(2) << std::setw(7) << v;
        }
        std::cout << std::endl;
    }
}

std::shared_ptr<MatrixCRS> MatrixCOO::toCRS() {
    std::stable_sort(data.begin(), data.end(),[](Triple t1, Triple t2) {return t1.row < t2.row;});

    auto crs = std::make_shared<MatrixCRS>(nRows, nCols, data.size());
    didx lastRow = 0;
    for (auto &trp: data) {
        didx row = trp.row;
        didx col = trp.col;
        dreal val = trp.val;
        crs->add(row-lastRow, col, val);
        lastRow = row;
    }
    return crs;
}

std::vector<Triple> MatrixCOO::vectorize() {
    return data;
}

void MatrixCOO::devectorize(std::vector<Triple> &vector) {
    data = vector;
    sumuped = false;
    ordered = false;
}

size_t MatrixCOO::getNumRows() {
    return nRows;
}

size_t MatrixCOO::getNumCols() {
    return nCols;
}

std::vector<size_t> MatrixCOO::getDimensions() {
    return std::vector<size_t>({nRows, nCols});
}

void MatrixCOO::setDimensions(std::vector<size_t> &&dimensions) {
    nRows = dimensions[0];
    nCols = dimensions[1];
}
