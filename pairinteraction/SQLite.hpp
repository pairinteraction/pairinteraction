/*
 * Copyright (c) 2016 Sebastian Weber, Henri Menke. All rights reserved.
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

#ifndef SQLITE_BINDINGS
#define SQLITE_BINDINGS

#include <iostream>
#include <memory>
#include <random>
#include <sstream>
#include <string>
#include <thread>
#include <type_traits>

#include <sqlite3.h>
#include <utility>

#include "utils.hpp"

namespace sqlite {

/** \brief SQLite error
 *
 * Exception to identify SQLite errors.
 * \code
 * try
 * {
 *   // code with sqlite
 * }
 * catch (sqlite::error &e)
 * {
 *   // handle exception
 * }
 * \endcode
 */
struct error : public std::exception {
public:
    /** \brief Constructor
     *
     * \param[in] err   error code
     * \param[in] msg   error message
     */
    explicit error(int err, std::string const &msg)
        : m_msg(std::string("SQLite error ") + std::to_string(err) + ": " + msg) {}

    /** \brief %what() function
     *
     * \returns error message
     */
    const char *what() const noexcept override { return m_msg.c_str(); }

private:
    std::string m_msg;
};

/** \brief SQLite Statement handle
 *
 * This object wraps an SQLite statement.  It is iterable in a while
 * loop with the step() function or in a range-based for loop.
 */
class statement {
    sqlite3 *m_db;
    std::unique_ptr<sqlite3_stmt, int (*)(sqlite3_stmt *)> m_stmt;
    std::string m_sql;
    bool m_prepared;
    bool m_valid;

    void handle_error(int err) {
        if (err != 0) {
            throw error(err, sqlite3_errstr(err));
        }
    }

public:
    /** \brief Constructor
     *
     * This is called in the query function of the the handle object.
     *
     * \param[in] db    unmanaged raw pointer to the sqlite database
     */
    explicit statement(sqlite3 *db) : statement{db, {}} {}

    /** \brief Constructor
     *
     * This is called in the query function of the the handle object.
     * This variant also sets a default query.
     *
     * \param[in] db    unmanaged raw pointer to the sqlite database
     * \param[in] sql   an initial query as a string
     */
    explicit statement(sqlite3 *db, std::string sql)
        : m_db{db}, m_stmt{nullptr, sqlite3_finalize}, m_sql{std::move(sql)},
          m_prepared{false}, m_valid{true} {}

    /** \brief Set the query string
     *
     * \param[in] sql   the query string
     * \throws sqlite::error
     */
    void set(std::string const &sql) {
        m_sql = sql;
        m_prepared = false;
    }

    /** \overload void set(std::string const &sql) */
    void set(std::stringstream const &ss) { return set(ss.str()); }

    /** \brief Prepare the statement
     *
     * Prepare the statement for stepping.
     *
     * \throws sqlite::error
     */
    void prepare() {
        sqlite3_stmt *pStmt; // is managed below
        auto err = sqlite3_prepare_v2(m_db, m_sql.c_str(), -1, &pStmt, nullptr);
        m_stmt.reset(pStmt);

        handle_error(err);

        m_prepared = true;
    }

    /** \brief Step the statement
     *
     * This steps the statement.
     *
     * \returns true if there is a row, false otherwise
     *
     * \throws sqlite::error
     */
    bool step() {
        if (!m_prepared) {
            handle_error(SQLITE_MISUSE);
        }
        if (!m_valid) {
            handle_error(SQLITE_DONE);
        }

        auto err = sqlite3_step(m_stmt.get());

        switch (err) {
        case SQLITE_DONE:
            m_valid = false;
            return false;
        case SQLITE_ROW:
            m_valid = true;
            return true;
        }

        m_valid = false;
        handle_error(err);
        return false;
    }

    /** \brief Reset the statement
     *
     * Reset the statement so you can rebind values.
     *
     * \throws sqlite::error
     */
    void reset() {
        handle_error(sqlite3_reset(m_stmt.get()));
        m_valid = true;
    }

    /** \brief Execute SQLite statements
     *
     * The SQL statements are passed to this function as a string and
     * are executed in-place.
     *
     * \param[in] sql   SQL statements
     * \throws sqlite::error
     */
    void exec(std::string const &sql) {
        set(sql);
        auto err = sqlite3_exec(m_db, m_sql.c_str(), nullptr, nullptr, nullptr);
        handle_error(err);
    }

    /** \brief Bind a value to a position in the query
     *
     * \param[in] where   position in the query (one-based)
     * \param[in] s       string to bind
     * \throws sqlite::error
     */
    void bind(int where, std::string const &s) {
        handle_error(sqlite3_bind_text(m_stmt.get(), where, s.c_str(), s.length(), SQLITE_STATIC));
    }

    /** \overload void bind(int where, std::string const &s) */
    void bind(int where, std::string &&s) {
        handle_error(
            sqlite3_bind_text(m_stmt.get(), where, s.c_str(), s.length(), SQLITE_TRANSIENT));
    }

    /** \overload void bind(int where, std::string const &s) */
    void bind(int where, double d) { handle_error(sqlite3_bind_double(m_stmt.get(), where, d)); }

    /** \overload void bind(int where, std::string const &s) */
    void bind(int where, int i) { handle_error(sqlite3_bind_int(m_stmt.get(), where, i)); }

#ifdef SCANNED_BY_DOXYGEN
    /** \brief Get value of a field
     *
     * Retrieve the value stored at position \p field as datatype \p
     * ReturnType.
     *
     * \tparam ReturnType   requested return type
     * \param[in]           field
     * \returns             the value stored in \p field as datatype \p
     * ReturnType
     */
    template <typename ReturnType>
    ReturnType get(int field);
#else
    template <typename T,
              typename = typename std::enable_if<std::is_floating_point<T>::value>::type>
    double get(int field) {
        return sqlite3_column_double(m_stmt.get(), field);
    }

    template <typename T, typename = typename std::enable_if<std::is_same<T, int>::value>::type>
    int get(int field) {
        return sqlite3_column_int(m_stmt.get(), field);
    }

    template <typename T,
              typename = typename std::enable_if<std::is_same<T, std::string>::value>::type>
    std::string get(int field) {
        return std::string(
            reinterpret_cast<char const *>(sqlite3_column_text(m_stmt.get(), field)));
    }
#endif

    /** \brief Statement iterator
     *
     * The iterator builds upon Boost Iterator Facade to implement an
     * iterator over the statement with minimal boilerplate.  We refer
     * to the official documentation for the internal functions:
     * http://www.boost.org/libs/iterator/
     */
    class iterator {
    private:
        statement *m_stmt;
        bool m_done;
        bool m_end;

    public:
        iterator() : m_stmt{nullptr}, m_done{true}, m_end{true} {}

        explicit iterator(statement *p) : m_stmt{p}, m_done{false}, m_end{false} { m_stmt->step(); }

        statement &operator*() const {
#ifndef NDEBUG
            if (m_end)
                throw std::out_of_range("iterator out of range");
#endif
            return *m_stmt;
        }
        iterator &operator++() {
            if (!m_done) {
                m_done = m_stmt->step();
            } else {
                m_end = true;
            }
            return *this;
        }
        bool operator!=(iterator const & /* unused */) const { return !m_end; }
    };

    /** \brief returns an iterator to the beginning
     *
     * With this, a statement can be iterated in a C++11 range-based for loop.
     *
     * \warning This is *not* zero overhead!  If you care about
     * performance use a while loop using step() instead.
     * \code
     * while(stmt.step())
     *     stmt.get<...>(...);
     * \endcode
     *
     * \returns iterator to the beginning
     */
    iterator begin() { return iterator{this}; }

    /** \brief returns an empty iterator
     *
     * The iterator for the SQLite statement is a forward iterator
     * which is implemented in terms of the step() function.  Since
     * the step() function will determine the end by itself, this
     * iterator is merely a sentinel.
     *
     * \returns empty iterator
     */
    iterator end() { return iterator{}; }
};

/** \brief SQLite Database handle
 *
 * This object handles the connection to the underlying SQLite databse.  It
 * provides a high-level object oriented RAII interface to the C-bindings of
 * SQLite3.
 */
class handle final {
    std::unique_ptr<sqlite3, decltype(&sqlite3_close)> m_db;
    int m_threshold;

    static int busy_handler(void *self, int num_prior_calls) {
        int thresh = static_cast<handle *>(self)->m_threshold;
        // Sleep if handler has been called less than num_prior_calls
        if (num_prior_calls < thresh) {
            std::this_thread::sleep_for(std::chrono::microseconds(utils::randint(2000, 20000)));
            return 1;
        }

        return 0; // Make sqlite3 return SQLITE_BUSY
    }

    void install_busy_handler() {
        auto err = sqlite3_busy_handler(m_db.get(), busy_handler, this);

        if (err != 0) {
            throw error(err, sqlite3_errmsg(m_db.get()));
        }
    }

public:
    /** \brief Conversion operator
     *
     * The handle object implicitly behaves like the underlying
     * pointer.
     *
     * \returns the raw database pointer
     */
    operator sqlite3 *() { return m_db.get(); }

    /** \brief Constructor
     *
     * It takes a filename and a bitmask of parameters to open the
     * database connection.  This is useful if you have to open the
     * connection, e.g., read-only.  By default it uses the standard
     * flags for opening a database connection.  See
     * https://www.sqlite.org/c3ref/open.html for details.
     *
     * \param[in] filename    fully-qualified filename of the database
     * \param[in] flags       options to open the database (optional)
     * \throws sqlite::error
     */
    explicit handle(std::string const &filename,
                    int flags = SQLITE_OPEN_READWRITE | SQLITE_OPEN_CREATE)
        : m_db{nullptr, sqlite3_close}, m_threshold{100000} {
        sqlite3 *tmp_db;
        auto err = sqlite3_open_v2(filename.c_str(), &tmp_db, flags, nullptr);
        m_db.reset(tmp_db);

        if (err != 0) {
            throw error(err, sqlite3_errmsg(m_db.get()));
        }

        install_busy_handler();
    }

    /** \brief Move constructor
     *
     * Reinstalls the busy handler.
     */
    handle(handle &&other) noexcept : m_db{std::move(other.m_db)}, m_threshold{other.m_threshold} {
        install_busy_handler();
    }

    /** \brief Move assignment operator
     *
     * Reinstalls the busy handler.
     */
    handle &operator=(handle &&other) noexcept {
        if (this != std::addressof(other)) // skip self-assignment
        {
            m_db = std::move(other.m_db);
            m_threshold = other.m_threshold;
            install_busy_handler();
        }
        return *this;
    }
};

} // namespace sqlite

#endif // SQLITE_BINDINGS
