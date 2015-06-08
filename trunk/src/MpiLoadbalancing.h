#ifndef MPILOADBALANCING_H
#define MPILOADBALANCING_H

#include "dtypes.h"

#include "MpiEnvironment.h"
#include "MpiVariables.h"
#include "Serializable.h"
#include "Vectorizable.h"

#include <math.h>
#include <vector>
#include <memory>

enum {WORKTAG = 0, DIETAG = 1};

class MpiLoadbalancing {
public:
    MpiLoadbalancing(std::shared_ptr<MpiEnvironment> mpi, int num);

protected:
    std::shared_ptr<MpiVariables> mympi;
    std::shared_ptr<MpiEnvironment> mpi;

    void runMaster(std::vector<std::shared_ptr<Serializable>> &dataIn, std::vector<std::shared_ptr<Serializable> > &dataOut);
    void runSlave(std::shared_ptr<Serializable> bufferSerializable, std::shared_ptr<Vectorizable> bufferVectorizable);

    virtual std::shared_ptr<Vectorizable> doMainprocessing(std::shared_ptr<Serializable> work) = 0; // all slaves within slave group
    virtual std::shared_ptr<Serializable> doPostprocessing(std::shared_ptr<Vectorizable> collection, std::shared_ptr<Serializable> work) = 0; // zeroth slave of slave group
    void run(std::vector<std::shared_ptr<Serializable> > &dataIn, std::vector<std::shared_ptr<Serializable>> &dataOut, std::shared_ptr<Serializable> bufferSerializable, std::shared_ptr<Vectorizable> bufferVectorizable);

private:
    int numSlaves;
    std::vector<int> vecSlave0;
    MPI_Datatype MPI_COSTUME;  
};

#endif
