#include "MpiEnvironment.h"

MpiEnvironment::MpiEnvironment(int argc, char **argv) : MpiVariables(Init(argc, argv)) {
}

MPI_Comm MpiEnvironment::Init(int argc, char **argv) {
    MPI_Init (&argc, &argv);
    return MPI_COMM_WORLD;
}

MpiEnvironment::~MpiEnvironment() {
    MPI_Finalize();
}
