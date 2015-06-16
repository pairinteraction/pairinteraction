#ifndef DTYPES_H
#define DTYPES_H

#include <vector>

#define REAL_T MPI_FLOAT
#define IDX_T MPI_UNSIGNED
#define BYTE_T MPI_UNSIGNED_CHAR

typedef float real_t;
typedef unsigned int idx_t;
typedef unsigned char byte_t;
typedef std::vector<byte_t> bytes_t;
typedef std::nullptr_t invalid_t;

struct Triple {
    Triple() : row(0), col(0), val(0) { }
    Triple(idx_t row, idx_t col, real_t val) : row(row), col(col), val(val) { }
    idx_t row;
    idx_t col;
    real_t val;
};

#endif
