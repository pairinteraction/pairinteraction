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

#ifndef MPIENVIRONMENT_H
#define MPIENVIRONMENT_H

#define OMPI_SKIP_MPICXX

#include "MpiVariables.h"
#include <inttypes.h>
#include <mpi.h>
#include <cassert>

class MpiEnvironment : public MpiVariables {
public:
    MpiEnvironment(int argc, char **argv);
    ~MpiEnvironment();
private:
    MPI_Comm Init(int argc, char **argv);
};

#endif
