#include "MpiEnvironment.h"

MpiEnvironment::MpiEnvironment(int argc, char **argv) : MpiVariables(Init(argc, argv)) {
}

const MPI::Intracomm& MpiEnvironment::Init(int argc, char **argv) {
    MPI::Init(argc, argv);
    assert(MPI::Is_initialized() == true);
    return MPI::COMM_WORLD;
}

MpiEnvironment::~MpiEnvironment() {
    MPI::Finalize();
    assert(MPI::Is_finalized() == true);
}
