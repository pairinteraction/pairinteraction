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

#ifndef SYSTEMONE_H
#define SYSTEMONE_H

#include "State.h"
#include "SystemBase.h"

#include <array>
#include <unordered_map>
#include <boost/functional/hash.hpp>

#ifndef SWIG
namespace std {
    template <>
    struct hash<std::array<int, 2>>
    {
        size_t operator()(const std::array<int, 2>& a) const
        {
            return boost::hash_value(a);
        }
    };
}
#endif

class SystemOne : public SystemBase<StateOne> {
public:
    SystemOne(std::string const& element, std::string cachedir);
    SystemOne(std::string const& element, std::string cachedir, bool memory_saving);
    SystemOne(std::string const& element);
    SystemOne(std::string const& element, bool memory_saving);

    const std::string& getElement() const;
    void setEfield(std::array<double, 3> field);
    void setBfield(std::array<double, 3> field);
    void setDiamagnetism(bool enable);

protected:
    void initializeBasis() override;
    void initializeInteraction() override;
    void addInteraction() override;
    void transformInteraction(const eigen_sparse_t &transformator) override;
    void deleteInteraction() override;

private:
    std::array<double, 3> efield, bfield;
    bool diamagnetism;
    std::string element;

    std::unordered_map<int, eigen_sparse_t> interaction_efield;
    std::unordered_map<int, eigen_sparse_t> interaction_bfield;
    std::unordered_map<std::array<int, 2>, eigen_sparse_t> interaction_diamagnetism;

    void changeToSphericalbasis(std::array<double, 3> field, std::unordered_map<int, double>& field_spherical);
    void changeToSphericalbasis(std::array<double, 3> field, std::unordered_map<int, std::complex<double>>& field_spherical);
};

#endif
