#ifndef MPILOADBALANCINGSIMPLE_H
#define MPILOADBALANCINGSIMPLE_H

#include "dtypes.h"

#include "MpiLoadbalancing.h"

#include <vector>
#include <memory>

class MpiLoadbalancingSimple : public MpiLoadbalancing  {
public:
    MpiLoadbalancingSimple(std::shared_ptr<MpiEnvironment> mpi, int num);

protected:
    virtual std::shared_ptr<Serializable> doProcessing(std::shared_ptr<Serializable> work) = 0;
    void run(std::vector<std::shared_ptr<Serializable> > &dataIn, std::vector<std::shared_ptr<Serializable>> &dataOut, std::shared_ptr<Serializable> bufferSerializable);

private:
    std::shared_ptr<Vectorizable> doMainprocessing(std::shared_ptr<Serializable> work);
    std::shared_ptr<Serializable> doPostprocessing(std::shared_ptr<Vectorizable> collection, std::shared_ptr<Serializable> work);
};

#endif
