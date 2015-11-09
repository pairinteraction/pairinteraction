#include "MpiVariables.h"

MpiVariables::MpiVariables(MPI::Intracomm world_) : world_(world_) {
    rank_ = world_.Get_rank();
    size_ = world_.Get_size();
}

int MpiVariables::rank() const {
    return rank_;
}

int MpiVariables::size() const {
    return size_;
}

const MPI::Intracomm& MpiVariables::world() const {
    return world_;
}
