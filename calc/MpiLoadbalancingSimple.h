/*
 * Copyright (c) 2016 Sebastian Weber, Henri Menke. All rights reserved.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

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
