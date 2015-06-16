#include "MatrixCRS.h"
#include "MatrixCOO.h"

MatrixCRS::MatrixCRS(size_t nRows, size_t nCols, size_t size) : nRows(nRows), nCols(nCols), ordered(true), sumuped(true) { // empty
    ptr.reserve(nRows+1);
    col.reserve(size);
    val.reserve(size);
    ptr.push_back(0);
}

MatrixCRS::MatrixCRS(bytes_t &bytes) { // from bytes
    deserialize(bytes);
}

MatrixCRS::MatrixCRS() {
}

void MatrixCRS::add(idx_t rIncrement, idx_t c, real_t v) {
    for (idx_t i = 0; i < rIncrement; ++i) {
        ptr.push_back(col.size());
    }
    col.push_back(c);
    val.push_back(v);
    ordered = false;
    sumuped = false;
}

void MatrixCRS::multiplyScalar(real_t &&scalar) {
    for (auto &v: val) {
        v *= scalar;
    }
}

void MatrixCRS::order() {
    if (!ordered) {
        while(ptr.size() < nRows+1) {
            ptr.push_back(col.size());
        }

        // Sort everything
        std::vector<size_t> indices;
        indices.reserve(col.size());
        for (size_t i = 0; i < col.size(); ++i) {
            indices.push_back(i);
        }

        for (idx_t r = 0; r < nRows; ++r) {
            idx_t ptrStart = (r < ptr.size()) ? ptr[r] : col.size();
            idx_t ptrEnd = (r+1 < ptr.size()) ? ptr[r+1] : col.size();
            std::sort(indices.begin()+ptrStart, indices.begin()+ptrEnd, [this](size_t i1, size_t i2) {return col[i1] < col[i2];});
        }

        auto colTmp(col);
        auto valTmp(val);
        for (size_t i = 0; i < indices.size(); ++i) {
            col[i] = colTmp[indices[i]];
            val[i] = valTmp[indices[i]];

        }

        ordered = true;
    }
}

void MatrixCRS::sumup() {
    order();

    if (!sumuped) {

        // Copy
        auto ptrTmp(ptr);
        auto colTmp(col);
        auto valTmp(val);
        ptr.clear();
        col.clear();
        val.clear();
        ptr.reserve(ptrTmp.size());
        col.reserve(colTmp.size());
        val.reserve(valTmp.size());
        ptr.push_back(0);

        // Sum doubled entries
        for (idx_t r = 0; r < nRows; ++r) {
            idx_t oldCol = 0;
            real_t sumVal = 0;

            for (idx_t p = ((r < ptrTmp.size()) ? ptrTmp[r] : colTmp.size()); p < ((r+1 < ptrTmp.size()) ? ptrTmp[r+1] : colTmp.size()); ++p) {
                idx_t c = colTmp[p];
                real_t v = valTmp[p];
                if (oldCol == c) {
                    sumVal += v;
                } else {
                    if (std::abs(sumVal) > 1e-12) {
                        col.push_back(oldCol);
                        val.push_back(sumVal);
                    }
                    sumVal = v;
                    oldCol = c;
                }
            }

            if (std::abs(sumVal) > 1e-12) {
                col.push_back(oldCol);
                val.push_back(sumVal);
            }

            ptr.push_back(col.size());
        }

        sumuped = true;
    }
}

void MatrixCRS::print() {
    order();
    sumup();

    std::cout << std::setiosflags(std::ios::right) << std::setiosflags(std::ios::fixed);
    for (idx_t r = 0; r < nRows; ++r) {
        idx_t p = (r < ptr.size()) ? ptr[r] : col.size();
        idx_t cc = col[p];

        for (idx_t c = 0; c < nCols; ++c) {
            real_t v = 0;
            if (p < ((r+1 < ptr.size()) ? ptr[r+1] : col.size()) && cc == c) {
                v = val[p];
                cc = col[++p];
            }
            std::cout << std::setprecision(2) << std::setw(7) << v;
        }
        std::cout << std::endl;
    }
}

std::shared_ptr<MatrixCOO>  MatrixCRS::toCOO() {
    auto coo = std::make_shared<MatrixCOO>(nRows, nCols, col.size());
    for (idx_t r = 0; r < nRows; ++r) {
        for (idx_t p = ((r < ptr.size()) ? ptr[r] : col.size()); p < ((r+1 < ptr.size()) ? ptr[r+1] : col.size()); ++p) {
            idx_t c = col[p];
            real_t v = val[p];
            coo->add(r,c,v);
        }
    }
    return coo;
}

bytes_t MatrixCRS::serialize() {
    order();

    size_t sizeHeader = sizeof(nRows)+sizeof(nCols)+sizeof(ordered)+sizeof(sumuped);
    size_t sizePtr = 3*sizeof(size_t)+sizeof(ptr[0])*ptr.size();
    size_t sizeCol = 3*sizeof(size_t)+sizeof(col[0])*col.size();
    size_t sizeVal = 3*sizeof(size_t)+sizeof(val[0])*val.size();
    bytes_t bytes(sizeHeader+sizePtr+sizeCol+sizeVal);

    auto pbytes = bytes.begin();
    serializeItem(pbytes,nRows);
    serializeItem(pbytes,nCols);
    serializeItem(pbytes,ordered);
    serializeItem(pbytes,sumuped);
    serializeItem(pbytes,ptr);
    serializeItem(pbytes,col);
    serializeItem(pbytes,val);

    return bytes;
}

void MatrixCRS::deserialize(bytes_t &bytes) {
    auto pbytes = bytes.begin();
    deserializeItem(pbytes,nRows);
    deserializeItem(pbytes,nCols);
    deserializeItem(pbytes,ordered);
    deserializeItem(pbytes,sumuped);
    deserializeItem(pbytes,ptr);
    deserializeItem(pbytes,col);
    deserializeItem(pbytes,val);
}

size_t MatrixCRS::getNumRows() {
    return nRows;
}

size_t MatrixCRS::getNumCols() {
    return nCols;
}

std::vector<size_t> MatrixCRS::getDimensions() {
    return std::vector<size_t>({nRows, nCols});
}

void MatrixCRS::setDimensions(std::vector<size_t> &&dimensions) {
    nRows = dimensions[0];
    nCols = dimensions[1];
}
