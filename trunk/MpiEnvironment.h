#ifndef MPIENVIRONMENT_H
#define MPIENVIRONMENT_H

#define OMPI_SKIP_MPICXX

#include "MpiVariables.h"
#include <inttypes.h>
#include <mpi.h>
#include <cassert>

class MpiEnvironment : public MpiVariables {
public:
    MpiEnvironment(int argc, char **argv);
    ~MpiEnvironment();
private:
    MPI_Comm Init(int argc, char **argv);
};

#endif
