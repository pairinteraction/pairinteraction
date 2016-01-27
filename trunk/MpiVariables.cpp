#include "MpiVariables.h"

MpiVariables::MpiVariables(MPI_Comm world_) : world_(world_) {
    MPI_Comm_rank(world_, &rank_);
    MPI_Comm_size(world_, &size_);
}

int MpiVariables::rank() const {
    return rank_;
}

int MpiVariables::size() const {
    return size_;
}

const MPI_Comm& MpiVariables::world() const {
    return world_;
}
