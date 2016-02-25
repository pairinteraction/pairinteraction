#ifndef MPIVARIABLES_H
#define MPIVARIABLES_H

#define OMPI_SKIP_MPICXX

#include <inttypes.h>
#include <mpi.h>

class MpiVariables {
public:
    MpiVariables(MPI_Comm world_);
    int rank() const;
    int size() const;
    const MPI_Comm& world() const;

private:
    int rank_;
    int size_;
    MPI_Comm world_;
};

#endif
