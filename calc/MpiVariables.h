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

#ifndef MPIVARIABLES_H
#define MPIVARIABLES_H

#define OMPI_SKIP_MPICXX

#include <inttypes.h>
#include <mpi.h>

class MpiVariables {
public:
    MpiVariables(MPI_Comm world_);
    int rank() const;
    int size() const;
    const MPI_Comm& world() const;

private:
    int rank_;
    int size_;
    MPI_Comm world_;
};

#endif
