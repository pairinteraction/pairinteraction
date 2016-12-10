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

#include "ConfParser.h"
#include <string>
#include <fstream>
#include <iostream>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

Configuration::value::value() : val() {
}
Configuration::value::value(std::stringstream v) {
    val.str(v.str());
}
std::string Configuration::value::str() const {
    return val.str();
}
void Configuration::value::str(std::string s) {
    val.str(s);
}
Configuration::value& Configuration::value::operator=(const Configuration::value& rhs) {
    val.str(std::string());
    val << rhs.str();
    return *this;
}
Configuration::iterator::entry::entry(const std::string key, Configuration::value& value) : key(key), value(value) {
}
Configuration::iterator::iterator(std::map<std::string, Configuration::value>::iterator itr) : itr(itr) {
}
bool Configuration::iterator::operator== (const iterator& other) const {
    return itr == other.itr;
}
bool Configuration::iterator::operator!= (const iterator& other) const {
    return itr != other.itr;
}
const Configuration::iterator& Configuration::iterator::operator++ () {
    ++itr;
    return *this;
}
Configuration::iterator::entry Configuration::iterator::operator* () {
    return iterator::entry(itr->first, itr->second);
}
Configuration::const_iterator::entry::entry(const std::string key, const Configuration::value& value) : key(key), value(value) {
}
Configuration::const_iterator::const_iterator(std::map<std::string, Configuration::value>::const_iterator itr) : itr(itr) {
}
bool Configuration::const_iterator::operator== (const const_iterator& other) const {
    return itr == other.itr;
}
bool Configuration::const_iterator::operator!= (const const_iterator& other) const {
    return itr != other.itr;
}
const Configuration::const_iterator& Configuration::const_iterator::operator++ () {
    ++itr;
    return *this;
}
Configuration::const_iterator::entry Configuration::const_iterator::operator* () {
    return const_iterator::entry(itr->first, itr->second);
}
int Configuration::load_from_json(std::string filename)
{
    using boost::property_tree::ptree;
    using boost::property_tree::json_parser::read_json;

    ptree pt;
    read_json(filename, pt);

    for (auto itr: pt)
    {
        params[itr.first] << itr.second.data();
    }

    return 0;
}
int Configuration::save_to_json(std::string filename)
{
    using boost::property_tree::ptree;
    using boost::property_tree::json_parser::write_json;

    ptree pt;

    for (auto p: *this) {
        pt.put(p.key,p.value.str());
    }

    write_json(filename, pt);

    return 0;
}
Configuration::value& Configuration::operator [](const std::string& key) {
    return params[key];
}
const Configuration::value& Configuration::operator [](const std::string& key) const{
    return params.at(key);
}
size_t Configuration::size() const {
    return params.size();
}
Configuration::iterator Configuration::find (const std::string& key) {
    return Configuration::iterator(params.find(key));
}

Configuration::const_iterator Configuration::find (const std::string& key) const {
    return Configuration::const_iterator(params.find(key));
}
Configuration::iterator Configuration::begin() {
    return Configuration::iterator(params.begin());
}
Configuration::iterator Configuration::end() {
    return Configuration::iterator(params.end());
}
Configuration::const_iterator Configuration::begin() const {
    return Configuration::const_iterator(params.begin());
}
Configuration::const_iterator Configuration::end() const {
    return Configuration::const_iterator(params.end());
}
bool Configuration::operator==(const Configuration &rhs) const {
    if (this->size() != rhs.size()) {
        return false;
    }
    for (auto p: *this) {
        if (rhs.find(p.key) == rhs.end() || rhs[p.key].str() != p.value.str()) {
            return false;
        }
    }
    return true;
}
Configuration& Configuration::operator+=(const Configuration& rhs) {
    for (auto p: rhs) {
        (*this)[p.key] = p.value;
    }
    return *this;
}
size_t Configuration::count(const std::string& key) const {
    return params.count(key);
}
Configuration::Configuration() {
}
