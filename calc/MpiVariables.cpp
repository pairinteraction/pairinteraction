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
