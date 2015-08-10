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

void MatrixCOO::add(idx_t row, idx_t col, real_t val) {
    data.push_back(Triple(row,col,val));
    ordered = false;
    sumuped = false;
}


void MatrixCOO::multiplyScalar(real_t &&scalar) {
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
        idx_t oldRow = 0;
        idx_t oldCol = 0;
        real_t sumVal = 0;
        for (auto &trp: tripleTmp) {
            idx_t row = trp.row;
            idx_t col = trp.col;
            real_t val = trp.val;

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

        data.shrink_to_fit();

        sumuped = true;
    }
}

void MatrixCOO::deleteZerocols() {
    // Search for columns with at least one non-zero entry
    std::vector<bool> nonzero(nCols,false);
    for (auto &trp: data) {
        if (std::abs(trp.val) > 1e-12) {
            nonzero[trp.col] = true;
        }
    }

    size_t i = 0;
    std::vector<size_t> mapping(nCols,false);
    for (size_t idx = 0; idx < nCols; ++idx) {
        if (nonzero[idx]) mapping[idx] = i++;
    }
    nCols = i++;

    // Copy
    auto tripleTmp(data);
    data.clear();
    data.reserve(tripleTmp.size());

    // Delete columns consisting of zero entries only
    for (auto &trp: tripleTmp) {
        if (std::abs(trp.val) > 1e-12) {
            data.push_back(Triple(trp.row,mapping[trp.col],trp.val));
        }
    }

    data.shrink_to_fit();
}

void MatrixCOO::deleteZerorows() {
    // Search for columns with at least one non-zero entry
    std::vector<bool> nonzero(nRows,false);
    for (auto &trp: data) {
        if (std::abs(trp.val) > 1e-12) {
            nonzero[trp.row] = true;
        }
    }

    size_t i = 0;
    std::vector<size_t> mapping(nRows,false);
    for (size_t idx = 0; idx < nRows; ++idx) {
        if (nonzero[idx]) mapping[idx] = i++;
    }
    nRows = i++;

    // Copy
    auto tripleTmp(data);
    data.clear();
    data.reserve(tripleTmp.size());

    // Delete columns consisting of zero entries only
    for (auto &trp: tripleTmp) {
        if (std::abs(trp.val) > 1e-12) {
            data.push_back(Triple(mapping[trp.row],trp.col,trp.val));
        }
    }

    data.shrink_to_fit();
}

void MatrixCOO::print() {
    order();
    sumup();

    idx_t n = 0;

    idx_t row = data[n].row;
    idx_t col = data[n].col;
    real_t val = data[n].val;

    std::cout << std::setiosflags(std::ios::right) << std::setiosflags(std::ios::fixed);
    for (idx_t r = 0; r < nRows; ++r) {
        for (idx_t c = 0; c < nCols; ++c) {
            real_t v = 0;
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

MatrixCRS MatrixCOO::toCRS() {
    std::stable_sort(data.begin(), data.end(),[](Triple t1, Triple t2) {return t1.row < t2.row;});

    MatrixCRS crs(nRows, nCols, data.size());
    for (auto &trp: data) {
        idx_t row = trp.row;
        idx_t col = trp.col;
        real_t val = trp.val;
        crs.add(row, col, val);
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
