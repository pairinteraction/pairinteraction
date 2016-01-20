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
        value(const Configuration::value &obj) {
            *this = obj;
        }
        std::string str() const;
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
        const Configuration::value& operator>>(T& rhs) const {
            std::stringstream tmp;
            tmp << val.str();
            tmp >> rhs;
            return *this;
        }
        template<typename T>
        Configuration::value& operator=(const T& rhs) {
            val.str(std::string());
            val << rhs;
            return *this;
        }
        Configuration::value& operator=(const Configuration::value& rhs);
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
        bool operator== (const iterator& other) const;
        const iterator& operator++ ();
        iterator::entry operator* ();
    private:
        std::map<std::string, Configuration::value>::iterator itr;
        int pos, nColumn;
        int rc;
    };

    class const_iterator {
    public:
        class entry {
        public:
            entry(const std::string key, const Configuration::value& value);
            const std::string key;
            const Configuration::value& value;
        };
        const_iterator(std::map<std::string, Configuration::value>::const_iterator itr);
        bool operator!= (const const_iterator& other) const;
        bool operator== (const const_iterator& other) const;
        const const_iterator& operator++ ();
        const_iterator::entry operator* ();
    private:
        std::map<std::string, Configuration::value>::const_iterator itr;
        int pos, nColumn;
        int rc;
    };

    int load_from_json(std::string filename);
    int save_to_json(std::string filename);
    Configuration::value &operator[](const std::string& key);
    const Configuration::value &operator[](const std::string& key) const;
    size_t size() const;
    Configuration::iterator find (const std::string& key);
    Configuration::const_iterator find (const std::string& key) const;
    Configuration::iterator begin();
    Configuration::iterator end();
    Configuration::const_iterator begin() const;
    Configuration::const_iterator end() const;
    bool operator==(const Configuration &rhs) const;
    size_t count(const std::string& key) const;
    Configuration();

private:
    std::map<std::string, Configuration::value> params;
};


#endif // CONF_PARSER_HPP
