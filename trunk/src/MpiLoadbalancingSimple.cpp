#include "MpiLoadbalancingSimple.h"

MpiLoadbalancingSimple::MpiLoadbalancingSimple(std::shared_ptr<MpiEnvironment> mpi, int num) : MpiLoadbalancing(mpi, num) {
}

void MpiLoadbalancingSimple::run(std::vector<std::shared_ptr<Serializable>> &dataIn, std::vector<std::shared_ptr<Serializable>> &dataOut, std::shared_ptr<Serializable> bufferSerializable) {
    if (mpi->rank() == 0) {
        this->runMaster(dataIn, dataOut);
    } else {
        this->runSlave(bufferSerializable, NULL);
    }
}

std::shared_ptr<Vectorizable> MpiLoadbalancingSimple::doMainprocessing(std::shared_ptr<Serializable> work) {
    (void)work;
    return NULL;
}

std::shared_ptr<Serializable> MpiLoadbalancingSimple::doPostprocessing(std::shared_ptr<Vectorizable> collection, std::shared_ptr<Serializable> work) {
    (void)collection;
    return doProcessing(work);
}
