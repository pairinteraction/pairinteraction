/*
 * Copyright (c) 2017 Sebastian Weber, Henri Menke. All rights reserved.
 *
 * This file is part of the pairinteraction library.
 *
 * The pairinteraction library is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * The pairinteraction library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with the pairinteraction library. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef CONF_PARSER_H
#define CONF_PARSER_H

#include <map>
#include <sstream>
#include <string>

/** \brief %Configuration storage
 *
 * This class stores the system configuration.  In essence it is a
 * wrapper around a `std::map` which `std::string` key and
 * `Configuration::value` value.  The wrapped value class has some
 * advanced features which make it better than a pure `std::string`
 * value.
 */
class Configuration {
private:
    class value {
        std::string m_value;

    public:
        value() = default;

        std::string const str() const { return m_value; }

        bool operator==(value const &rhs) const { return this->m_value == rhs.m_value; }

        template <typename T>
        value &operator<<(T const &rhs) {
            std::ostringstream os;
            os << rhs;
            m_value = os.str();
            return *this;
        }

        value &operator<<(std::string const &rhs) {
            m_value = rhs;
            return *this;
        }

        value &operator<<(char const *rhs) {
            m_value = std::string{rhs};
            return *this;
        }

        value &operator>>(int &rhs) {
            rhs = std::stoi(m_value);
            return *this;
        }

        value &operator>>(size_t &rhs) {
            rhs = std::stoul(m_value);
            return *this;
        }

        value &operator>>(float &rhs) {
            rhs = std::stof(m_value);
            return *this;
        }

        value &operator>>(double &rhs) {
            rhs = std::stod(m_value);
            return *this;
        }

        value &operator>>(std::string &rhs) {
            rhs = m_value;
            return *this;
        }

        value const &operator>>(int &rhs) const {
            rhs = std::stoi(m_value);
            return *this;
        }

        value const &operator>>(size_t &rhs) const {
            rhs = std::stoul(m_value);
            return *this;
        }

        value const &operator>>(float &rhs) const {
            rhs = std::stof(m_value);
            return *this;
        }

        value const &operator>>(double &rhs) const {
            rhs = std::stod(m_value);
            return *this;
        }

        value const &operator>>(std::string &rhs) const {
            rhs = m_value;
            return *this;
        }

        value &operator<<(value const &rhs) {
            m_value = rhs.m_value;
            return *this;
        }

        value &operator>>(value &rhs) {
            rhs.m_value = m_value;
            return *this;
        }

        value const &operator>>(value &rhs) const {
            rhs.m_value = m_value;
            return *this;
        }
    };

public:
    /** \brief Load configuration from JSON file
     *
     * Loads the configuration from a JSON file using Boost's ptree.
     *
     * \param[in] filename    Path to the JSON file
     */
    void load_from_json(std::string const &filename);

    /** \brief Save configuration to JSON file
     *
     * Saves the configuration to a JSON file using Boost's ptree.
     *
     * \param[in] filename    Path to the JSON file
     */
    void save_to_json(std::string const &filename) const;

    /** \brief Constructor */
    Configuration() = default;

    /** \brief Number of elements matching specific key
     *
     * Returns the number of elements matching the key.  This can only
     * be zero or one, because the wrapped `std::map` does not allow
     * duplicates.
     *
     * \param[in] key    Key
     * \returns Number of elements machting key
     */
    size_t count(std::string const &key) const { return params.count(key); }

    /** \brief Number of elements
     *
     * \returns Number of elements
     */
    size_t size() const { return params.size(); }

    /** \brief Append operator
     *
     * This will append the right hand side to the left hand side.  If
     * a key exists on both sides, it will be overwritten with the
     * value from the right hand side.
     *
     * \returns Configuration appended to
     */
    Configuration &operator+=(Configuration const &rhs) {
        for (auto const &p : rhs) {
            (*this)[p.first] = p.second;
        }
        return *this;
    }

    /** \brief Element access operator
     *
     * This returns a reference to the element associated to the key
     * and will insert it if no such key exists.
     *
     * \returns Reference to the value
     */
    value &operator[](std::string const &key) { return params[key]; }

    /** \brief Constant element access operator
     *
     * This returns a copy of the element associated to the key.  It
     * throws if no such key exists.
     *
     * \returns Copy of the value
     * \throws std::out_of_range
     */
    value operator[](std::string const &key) const { return params.at(key); }

    /** \brief Comparison operator
     *
     * Compares two configurations.
     *
     * \returns Truth value for equality
     */
    bool operator==(Configuration const &rhs) const { return this->params == rhs.params; }

    /** \brief Iterator type
     *
     * Instead of implementing an own iterator we use the existing one
     * of the underlying `std::map`.
     */
    typedef std::map<std::string, value>::iterator iterator;

    /** \brief Constant iterator type
     *
     * Instead of implementing an own constant iterator we use the
     * existing one of the underlying `std::map`.
     */
    typedef std::map<std::string, value>::const_iterator const_iterator;

    /** \brief Iterator pointing to the beginning
     *
     * \returns iterator
     */
    iterator begin() { return params.begin(); }

    /** \brief Iterator pointing to the end
     *
     * \returns iterator
     */
    iterator end() { return params.end(); }

    /** \brief Constant iterator pointing to the beginning
     *
     * \returns constant iterator
     */
    const_iterator begin() const { return params.cbegin(); }

    /** \brief Constant iterator pointing to the end
     *
     * \returns constant iterator
     */
    const_iterator end() const { return params.cend(); }

    /** \brief Find an element associated to a key
     *
     * Attempts to find a key in the map and returns an iterator
     * pointing to the element.  If no such key exists it returns
     * end().
     *
     * \returns iterator
     */
    iterator find(std::string const &key) { return params.find(key); }

    /** \brief Constant version of find(std::string const&)
     *
     * \returns constant iterator
     */
    const_iterator find(std::string const &key) const { return params.find(key); }

private:
    std::map<std::string, value> params;
};

#endif // CONF_PARSER_H
