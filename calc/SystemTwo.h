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

#ifndef SYSTEMTWO_H
#define SYSTEMTWO_H

#include "State.h"
#include "SystemBase.h"
#include "SystemOne.h"

class SystemTwo : public SystemBase<StateTwo> {
public:
    SystemTwo(const SystemOne &b1, const SystemOne &b2);
    std::vector<StateOne> getStatesFirst();
    std::vector<StateOne> getStatesSecond();
protected:
    void initialize() override;
private:
    SystemOne basis1; // is needed in the initialize method and afterwards deleted
    SystemOne basis2; // is needed in the initialize method and afterwards deleted
};

#endif
