#ifndef DTYPES_H
#define DTYPES_H

#define DREAL MPI_FLOAT
#define DIDX MPI_UNSIGNED

typedef float dreal;
typedef unsigned int didx;

struct triple {
    triple() : row(0), col(0), val(0) { }
    triple(didx row, didx col, dreal val) : row(row), col(col), val(val) { }
    didx row;
    didx col;
    dreal val;
};

#endif
