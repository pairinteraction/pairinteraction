#ifndef DTYPES_H
#define DTYPES_H

#define DREAL MPI_FLOAT
#define DIDX MPI_UNSIGNED

typedef float dreal;
typedef unsigned int didx;

struct Triple {
    Triple() : row(0), col(0), val(0) { }
    Triple(didx row, didx col, dreal val) : row(row), col(col), val(val) { }
    didx row;
    didx col;
    dreal val;
};

#endif
