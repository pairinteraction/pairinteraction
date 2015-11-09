#ifndef MPIENVIRONMENT_H
#define MPIENVIRONMENT_H

#define _MPI_CPP_BINDINGS

#include "MpiVariables.h"
#include <mpi.h>
#include <cassert>

class MpiEnvironment : public MpiVariables {
public:
    MpiEnvironment(int argc, char **argv);
    ~MpiEnvironment();
private:
    const MPI::Intracomm& Init(int argc, char **argv);
};

#endif
