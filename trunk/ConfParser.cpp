#include "ConfParser.hpp"
#include <string>
#include <fstream>
#include <iostream>
#include "json/json.h"

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
    Json::Value datadict;
    Json::Reader reader;
    std::ifstream conffile(filename);
    if ( !reader.parse(conffile, datadict) )
    {
        std::cerr << reader.getFormattedErrorMessages();
        conffile.close();
        return 1;
    }
    conffile.close();

    for(Json::Value::iterator itr = datadict.begin(); itr !=datadict.end(); ++itr) {
        std::string s = (*itr).toStyledString();
        s.erase(remove( s.begin(), s.end(), '\"' ), s.end()); // TODO remove hack
         s.erase(remove( s.begin(), s.end(), '\n' ), s.end()); // TODO remove hack
        params[itr.key().asString()] << s;
    }

    return 0;
}
int Configuration::save_to_json(std::string filename)
{
    Json::StyledWriter styledWriter;
    Json::Value datadict;
    std::ofstream conffile(filename);
    for (auto p: *this) {
        datadict[p.key] = p.value.str();
    }
    conffile << styledWriter.write(datadict);
    conffile.close();

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
    return Configuration::iterator(std::move(params.find(key)));
}

Configuration::const_iterator Configuration::find (const std::string& key) const {
    return Configuration::const_iterator(std::move(params.find(key)));
}
Configuration::iterator Configuration::begin() {
    return Configuration::iterator(std::move(params.begin()));
}
Configuration::iterator Configuration::end() {
    return Configuration::iterator(std::move(params.end()));
}
Configuration::const_iterator Configuration::begin() const {
    return Configuration::const_iterator(std::move(params.begin()));
}
Configuration::const_iterator Configuration::end() const {
    return Configuration::const_iterator(std::move(params.end()));
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
