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

#ifndef CONF_PARSER_H
#define CONF_PARSER_H

#include <string>
#include <sstream>
#include <map>
#include "dtypes.h"

class Configuration {
public:
    class value
    {
        std::string m_value;
    public:
        value() : m_value() {};

        std::string const str() const
        {
            return m_value;
        }

        template < typename T >
        friend value& operator<<(value& v, T const& rhs)
        {
            std::stringstream ss;
            ss << rhs;
            ss >> v.m_value;
            return v;
        }
        
        template < typename T >
        friend value& operator>>(value& v, T& rhs)
        {
            std::stringstream ss;
            ss << v.m_value;
            ss >> rhs;
            return v;
        }

        template < typename T >
        friend value const& operator>>(value const& v, T& rhs)
        {
            std::stringstream ss;
            ss << v.m_value;
            ss >> rhs;
            return v;
        }

        friend std::stringstream& operator<<(std::stringstream& ss, value const& rhs)
        {
            ss << rhs.m_value;
            return ss;
        }

        friend std::stringstream& operator>>(std::stringstream& ss, value& rhs)
        {
            ss >> rhs.m_value;
            return ss;
        }
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

    void load_from_json(std::string filename);
    void save_to_json(std::string filename);
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
    Configuration& operator+=(const Configuration& rhs);
    size_t count(const std::string& key) const;
    Configuration();

private:
    std::map<std::string, Configuration::value> params;
};


#endif // CONF_PARSER_H
