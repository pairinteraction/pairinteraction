#ifndef CONF_PARSER_HPP
#define CONF_PARSER_HPP

#include <string>
#include <sstream>
#include <map>
#include "dtypes.h"

class Configuration {
public:
    class value {
    public:
        value();
        value(std::stringstream val);
        std::string str();
        void str(std::string s);
        template<typename T>
        Configuration::value& operator<<(const T& rhs) {
            val << rhs;
            return *this;
        }
        template<typename T>
        Configuration::value& operator>>(T& rhs) {
            val >> rhs;
            return *this;
        }
        template<typename T>
        Configuration::value& operator=(const T& rhs) {
            val.str(std::string());
            val << rhs;
            return *this;
        }
        Configuration::value& operator=(Configuration::value& rhs);
    private:
        std::stringstream val;
    };

    class iterator {
    public:
        class entry {
        public:
            entry(const std::string key, Configuration::value& value);
            const std::string key;
            Configuration::value& value;
        };
        iterator(std::map<std::string, Configuration::value>::iterator itr);
        bool operator!= (const iterator& other) const;
        const iterator& operator++ ();
        iterator::entry operator* ();
    private:
        std::map<std::string, Configuration::value>::iterator itr;
        int pos, nColumn;
        int rc;
    };

    int load_from_json(std::string filename);
    int save_to_json(std::string filename);
    Configuration::value &operator[](std::string key);
    size_t size() const;
    Configuration::iterator begin();
    Configuration::iterator end();
    bool operator==(Configuration &rhs);

private:
    std::map<std::string, Configuration::value> params;
};


#endif // CONF_PARSER_HPP
