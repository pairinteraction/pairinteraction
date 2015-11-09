#ifndef MPILOADBALANCINGSIMPLE_H
#define MPILOADBALANCINGSIMPLE_H

#include "dtypes.h"

#include "MpiLoadbalancing.h"

#include <vector>
#include <memory>

template <class TIn, class TOut>
class MpiLoadbalancingSimple : public MpiLoadbalancing<TIn, TOut>{
public:
    MpiLoadbalancingSimple(std::shared_ptr<MpiEnvironment> mpi, int num) : MpiLoadbalancing<TIn, TOut>(mpi, num) {
    }

protected:
    void runSlave() {
        if (this->mympi->rank() != 0) return;

        // === Do work ===
        while (true) {
            // --- Receive task from master ---
            Message message_in = this->ReceiveFromMaster();

            // --- Stop working in case of DIETAG ---
            if (message_in.size == -1) break;

            // --- Do calculations ---
            auto bufferSerializable = std::make_shared<TIn>();
            bufferSerializable->deserialize(message_in.data);
            auto results = doProcessing(std::move(bufferSerializable))->serialize();

            // --- Send results to master ---
            Message message_out(results, results.size(), 0);
            this->SendToMaster(message_out);
        }
    }

    virtual std::shared_ptr<TOut> doProcessing(std::shared_ptr<TIn> work) = 0;
};

#endif
