#ifndef VECTORIZABLE_H
#define VECTORIZABLE_H

#include "dtypes.h"

#include <vector>

class Vectorizable {
public:
    virtual std::vector<triple> vectorize() = 0;
    virtual void devectorize(std::vector<triple> &vector) = 0;
};

#endif
