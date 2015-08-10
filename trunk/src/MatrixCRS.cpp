#include "MatrixCRS.h"
#include "MatrixCOO.h"

MatrixCRS::MatrixCRS(size_t nRows, size_t nCols, size_t size) : nRows(nRows), nCols(nCols), r_last(0), ordered(true), sumuped(true) { // empty
    ptr_.reserve(nRows+1);
    col_.reserve(size);
    row_.reserve(size);
    val_.reserve(size);
    ptr_.push_back(0);
}

MatrixCRS::MatrixCRS(bytes_t &bytes) { // from bytes
    deserialize(bytes);
}

MatrixCRS::MatrixCRS() {
}

void MatrixCRS::init(size_t nRows, size_t nCols, size_t size) {
    this->nRows = nRows;
    this->nCols = nCols;
    this->ordered = true;
    this->sumuped = true;

    ptr_.clear();
    col_.clear();
    row_.clear();
    val_.clear();
    ptr_.clear();

    ptr_.reserve(nRows+1);
    col_.reserve(size);
    row_.reserve(size);
    val_.reserve(size);
    ptr_.push_back(0);
}

/*void MatrixCRS::add(idx_t rIncrement, idx_t c, real_t v) {
    for (idx_t i = 0; i < rIncrement; ++i) {
        ptr.push_back(col.size());
    }
    col.push_back(c);
    val.push_back(v);
    ordered = false;
    sumuped = false;
}*/

void MatrixCRS::add(idx_t r, idx_t c, real_t v) {
    idx_t rIncrement = r - r_last;
    r_last = r;

    for (idx_t i = 0; i < rIncrement; ++i) {
        ptr_.push_back(col_.size());
    }

    col_.push_back(c);
    row_.push_back(r);
    val_.push_back(v);
    ordered = false;
    sumuped = false;
}

void MatrixCRS::multiplyScalar(real_t &&scalar) {
    for (auto &v: val_) {
        v *= scalar;
    }
}

void MatrixCRS::order() {
    if (!ordered) {
        while(ptr_.size() < nRows+1) {
            ptr_.push_back(col_.size());
        }

        // Sort everything
        std::vector<size_t> indices;
        indices.reserve(col_.size());
        for (size_t i = 0; i < col_.size(); ++i) {
            indices.push_back(i);
        }

        for (idx_t r = 0; r < nRows; ++r) {
            idx_t ptrStart = (r < ptr_.size()) ? ptr_[r] : col_.size();
            idx_t ptrEnd = (r+1 < ptr_.size()) ? ptr_[r+1] : col_.size();
            std::sort(indices.begin()+ptrStart, indices.begin()+ptrEnd, [this](size_t i1, size_t i2) {return col_[i1] < col_[i2];});
        }

        auto colTmp(col_);
        auto valTmp(val_);
        for (size_t i = 0; i < indices.size(); ++i) {
            col_[i] = colTmp[indices[i]];
            val_[i] = valTmp[indices[i]];
        }

        buildRow();

        ordered = true;
    }
}

void MatrixCRS::sumup() {
    order();

    if (!sumuped) {

        // Copy
        auto ptrTmp(ptr_);
        auto colTmp(col_);
        auto valTmp(val_);
        ptr_.clear();
        col_.clear();
        val_.clear();
        ptr_.reserve(ptrTmp.size());
        col_.reserve(colTmp.size());
        val_.reserve(valTmp.size());
        ptr_.push_back(0);

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
                        col_.push_back(oldCol);
                        val_.push_back(sumVal);
                    }
                    sumVal = v;
                    oldCol = c;
                }
            }

            if (std::abs(sumVal) > 1e-12) {
                col_.push_back(oldCol);
                val_.push_back(sumVal);
            }

            ptr_.push_back(col_.size());
        }

        buildRow();

        sumuped = true;
    }
}

void MatrixCRS::print() {
    order();
    sumup();

    std::cout << std::setiosflags(std::ios::right) << std::setiosflags(std::ios::fixed);
    for (idx_t r = 0; r < nRows; ++r) {
        idx_t p = (r < ptr_.size()) ? ptr_[r] : col_.size();
        idx_t cc = col_[p];

        for (idx_t c = 0; c < nCols; ++c) {
            real_t v = 0;
            if (p < ((r+1 < ptr_.size()) ? ptr_[r+1] : col_.size()) && cc == c) {
                v = val_[p];
                cc = col_[++p];
            }
            std::cout << std::setprecision(2) << std::setw(7) << v;
        }
        std::cout << std::endl;
    }
}

RowCRS MatrixCRS::row(size_t row) {
    idx_t start = ((row < ptr_.size()) ? ptr_[row] : col_.size());
    idx_t stop = ((row+1 < ptr_.size()) ? ptr_[row+1] : col_.size());

    std::vector<idx_t> cols;
    std::vector<real_t> vals;
    cols.reserve(stop-start);
    vals.reserve(stop-start);

    for (idx_t p = start; p < stop; ++p) {
        cols.push_back(col_[p]);
        vals.push_back(val_[p]);
    }
    return RowCRS(row, cols, vals);
}

MatrixCRS MatrixCRS::transformBackward(MatrixCRS &basis) {
    // --- transformation ---
    MatrixCRS transformed(basis.getNumRows(),basis.getNumRows(),basis.getNumRows()*basis.getNumRows());

    std::unordered_map<idx_t, real_t> colvals_temp;
    std::unordered_map<idx_t, real_t> colvals_transformed;

    for (idx_t row_basis1=0; row_basis1<basis.getNumRows();++row_basis1) {
        colvals_temp.clear();
        colvals_transformed.clear();

        // Temp_ij = Basis_ik * Original_kj for all j  (i->k->j)
        for (auto triple_basis1 : basis.row(row_basis1)) { //k-loop
            for (auto triple_original : this->row(triple_basis1.col)) { //j-loop
                real_t val = triple_basis1.val*triple_original.val;
                if (colvals_temp.count(triple_original.col)) colvals_temp[triple_original.col] += val;
                else colvals_temp[triple_original.col] = val;
            }
        }

        // Transformed_il = Temp_ij * Basis.T_jl = Temp_ij * Basis_lj for all l (l->j)
        for (idx_t row_basis2=0;row_basis2<basis.getNumRows();++row_basis2) { //l-loop
            for (auto triple_basis2 : basis.row(row_basis2)) { //j-loop
                if (colvals_temp.count(triple_basis2.col)) {
                    real_t val_temp = colvals_temp[triple_basis2.col];
                    real_t val = val_temp*triple_basis2.val;
                    if (colvals_transformed.count(row_basis2)) colvals_transformed[row_basis2] += val;
                    else colvals_transformed[row_basis2] = val;
                }
            }
        }

        // save the result
        for (auto &item_transformed: colvals_transformed) {
            if (std::abs(item_transformed.second) > 1e-12) {
                transformed.add(row_basis1, item_transformed.first, item_transformed.second);
            }
        }
    }

    // --- save space ---
    // data.shrink_to_fit();

    return transformed;
}

MatrixCRS MatrixCRS::transformForward(MatrixCRS &basisoriginal) {
    // --- transpose basisoriginal ---
    MatrixCRS basis(basisoriginal.getNumCols(),basisoriginal.getNumRows(),basisoriginal.size());

    std::vector<std::vector<Triple> > transposer(basis.getNumRows());
    for (auto triple : basisoriginal) {
        transposer[triple.col].push_back(Triple(triple.row,triple.col,triple.val));
    }

    for (idx_t row=0; row<basis.getNumRows();++row) {
        for (auto triple : transposer[row]) {
            basis.add(triple.col,triple.row,triple.val);
        }
    }

    // --- transformation ---
    return transformBackward(basis);
}

MatrixCOO MatrixCRS::toCOO() {
    MatrixCOO coo(nRows, nCols, col_.size());
    for (idx_t r = 0; r < nRows; ++r) {
        for (idx_t p = ((r < ptr_.size()) ? ptr_[r] : col_.size()); p < ((r+1 < ptr_.size()) ? ptr_[r+1] : col_.size()); ++p) {
            idx_t c = col_[p];
            real_t v = val_[p];
            coo.add(r,c,v);
        }
    }
    return coo;
}

bytes_t MatrixCRS::serialize() {
    order();

    size_t sizeHeader = sizeof(nRows)+sizeof(nCols)+sizeof(ordered)+sizeof(sumuped);
    size_t sizePtr = 3*sizeof(size_t)+sizeof(ptr_[0])*ptr_.size();
    size_t sizeCol = 3*sizeof(size_t)+sizeof(col_[0])*col_.size();
    size_t sizeVal = 3*sizeof(size_t)+sizeof(val_[0])*val_.size();
    bytes_t bytes(sizeHeader+sizePtr+sizeCol+sizeVal);

    auto pbytes = bytes.begin();
    serializeItem(pbytes,nRows);
    serializeItem(pbytes,nCols);
    serializeItem(pbytes,ordered);
    serializeItem(pbytes,sumuped);
    serializeItem(pbytes,ptr_);
    serializeItem(pbytes,col_);
    serializeItem(pbytes,val_);

    return bytes;
}

void MatrixCRS::deserialize(bytes_t &bytes) {
    auto pbytes = bytes.begin();
    deserializeItem(pbytes,nRows);
    deserializeItem(pbytes,nCols);
    deserializeItem(pbytes,ordered);
    deserializeItem(pbytes,sumuped);
    deserializeItem(pbytes,ptr_);
    deserializeItem(pbytes,col_);
    deserializeItem(pbytes,val_);

    buildRow();
}

size_t MatrixCRS::getNumRows() {
    return nRows;
}

size_t MatrixCRS::getNumCols() {
    return nCols;
}

size_t MatrixCRS::size() {
    return val_.size();
}

std::vector<size_t> MatrixCRS::getDimensions() {
    return std::vector<size_t>({nRows, nCols});
}

void MatrixCRS::setDimensions(std::vector<size_t> &&dimensions) {
    nRows = dimensions[0];
    nCols = dimensions[1];
}

void MatrixCRS::buildRow() {
    row_.clear();
    row_.reserve(col_.size());
    for (idx_t r = 0; r < nRows; ++r) {
        for (idx_t p = ((r < ptr_.size()) ? ptr_[r] : col_.size()); p < ((r+1 < ptr_.size()) ? ptr_[r+1] : col_.size()); ++p) {
            row_.push_back(r);
        }
    }
}
