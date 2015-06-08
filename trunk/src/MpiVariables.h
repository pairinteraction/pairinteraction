#ifndef MPIVARIABLES_H
#define MPIVARIABLES_H

#define _MPI_CPP_BINDINGS

#include <mpi.h>

class MpiVariables {
public:
    MpiVariables(MPI::Intracomm world_);
    int rank() const;
    int size() const;
    const MPI::Intracomm& world() const;

private:
    int rank_;
    int size_;
    MPI::Intracomm world_;
};

#endif
