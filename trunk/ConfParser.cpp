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
std::string Configuration::value::str() {
    return val.str();
}
void Configuration::value::str(std::string s) {
    val.str(s);
}
Configuration::value& Configuration::value::operator=(Configuration::value& rhs) {
    val.str(std::string());
    val << rhs.str();
    return *this;
}
Configuration::iterator::entry::entry(const std::string key, Configuration::value& value) : key(key), value(value) {
}
Configuration::iterator::iterator(std::map<std::string, Configuration::value>::iterator itr) : itr(itr) {
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
        params[itr.key().asString()] << itr->asString();
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
Configuration::value& Configuration::operator [](std::string key) {
    return params[key];
}
size_t Configuration::size() const {
    return params.size();
}
Configuration::iterator Configuration::begin() {
    return Configuration::iterator(std::move(params.begin()));
}
Configuration::iterator Configuration::end() {
    return Configuration::iterator(std::move(params.end()));
}
bool Configuration::operator==(Configuration &rhs) { // TODO bool Configuration::operator==(const Configuration &rhs) const {
    if (this->size() != rhs.size()) {
        return false;
    }
    for (auto p: *this) {
        if (rhs[p.key].str() != p.value.str()) {
            return false;
        }
    }
    return true;
}
