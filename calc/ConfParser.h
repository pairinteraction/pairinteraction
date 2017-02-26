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

    void load_from_json(std::string const& filename);
    void save_to_json(std::string const& filename);

    Configuration() {};

    size_t count(std::string const& key) const
    {
        return params.count(key);
    }

    size_t size() const
    {
        return params.size();
    }

    Configuration& operator+=(Configuration const& rhs)
    {
        for (auto const& p : rhs)
        {
            (*this)[p.first] = p.second;
        }
        return *this;
    }

    value& operator[](std::string const& key)
    {
        return params[key];
    }

    value operator[](std::string const& key) const
    {
        return params.at(key);
    }

    bool operator==(Configuration const& rhs) const
    {
        if (this->size() != rhs.size())
        {
            return false;
        }
        for (auto const& p : *this)
        {
            if ( rhs.find(p.first) == rhs.end() ||
                 rhs[p.first].str() != p.second.str() )
            {
                return false;
            }
        }
        return true;
    }

    typedef std::map < std::string, value > ::iterator iterator;
    typedef std::map < std::string, value > ::const_iterator const_iterator;

    iterator begin()
    {
        return params.begin();
    }

    iterator end()
    {
        return params.end();
    }

    const_iterator begin() const
    {
        return params.cbegin();
    }

    const_iterator end() const
    {
        return params.cend();
    }

    iterator find(std::string const& key)
    {
        return params.find(key);
    }

    const_iterator find(std::string const& key) const
    {
        return params.find(key);
    }

private:
    std::map<std::string, Configuration::value> params;
};


#endif // CONF_PARSER_H
