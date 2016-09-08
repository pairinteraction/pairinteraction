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

#ifndef MPILOADBALANCINGCOMPLEX_H
#define MPILOADBALANCINGCOMPLEX_H

#include "MpiLoadbalancing.h"

#include <vector>
#include <memory>


template <class TIn, class TOut, class TVec>
class MpiLoadbalancingComplex : public MpiLoadbalancing<TIn, TOut>{
public:
    MpiLoadbalancingComplex(std::shared_ptr<MpiEnvironment> mpi, int num) : MpiLoadbalancing<TIn, TOut>(mpi, num) {
    }

protected:
    void runSlave() {
        // === Do work ===
        while (true) {
            Message message_in;

            // --- Receive message from master ---
            if (this->mympi->rank() == 0) {
                message_in = this->ReceiveFromMaster();
            }

            // --- Receive message from slave0 ---
            if (this->mympi->size() > 1) {
                // Broadcast message size
                MPI_Bcast(&message_in.size, 1, MPI_INT, 0, this->mympi->world());

                if (message_in.size != -1) {
                    if (this->mympi->rank() != 0 ) {
                        // Allocate buffer
                        message_in.data.resize(message_in.size);
                    }

                    // Broadcast message
                    this->mympi->world().Bcast(&message_in.data[0], message_in.size, BYTE_T, 0);
                }
            }

            // --- Stop working in case of DIETAG ---
            if (message_in.size == -1) break;

            // --- Do calculations ---
            auto bufferSerializable = std::make_shared<TIn>();
            bufferSerializable->deserialize(message_in.data);
            auto resultunvectorized = doMainprocessing(bufferSerializable);

            // --- Send results to slave0 ---
            std::vector<Triple> collection;
            auto result = resultunvectorized->vectorize();

            if (this->mympi->size() > 1) {
                collection = this->Gathervector(result);
            } else {
                collection = result;
            }

            // --- Do postprocessing and send results to master ---
            if (this->mympi->rank() == 0) {
                // Do postprocessing
                auto bufferVectorizable = std::make_shared<TVec>();
                bufferVectorizable->devectorize(collection);
                auto results = doPostprocessing(std::move(bufferVectorizable), std::move(bufferSerializable))->serialize();

                // Send results to master
                Message message_out(results, results.size(), 0);
                this->SendToMaster(message_out);
            }
        }
    }

    virtual std::shared_ptr<TVec> doMainprocessing(std::shared_ptr<TIn> work) = 0; // all slaves within slave group
    virtual std::shared_ptr<TOut> doPostprocessing(std::shared_ptr<TVec> collection, std::shared_ptr<TIn> work) = 0; // zeroth slave of slave group
};

#endif
