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

#ifndef BASIS_H
#define BASIS_H

#include "State.h"
#include "BasisBase.h"

class BasisOne : public Basis<StateOne> {
public:
    BasisOne(std::string const& element);
    const std::string& getElement() const;
protected:
    void initialize() override;
private:
    std::string element;
};

class BasisTwo : public Basis<StateTwo> {
public:
    BasisTwo(const BasisOne &b1, const BasisOne &b2);
    BasisOne getFirstBasis() const;
    BasisOne getSecondBasis() const;
protected:
    void initialize() override;
private:
    BasisOne basis1; // is needed in the initialize method and afterwards deleted
    BasisOne basis2; // is needed in the initialize method and afterwards deleted
};

#endif
