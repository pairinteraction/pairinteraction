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

#ifndef SQLITE_BINDINGS
#define SQLITE_BINDINGS

#include <functional>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <type_traits>

#include <boost/iterator/iterator_facade.hpp>
#include <sqlite3.h>

namespace sqlite
{

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
        : m_msg(std::string("SQLite error ") + std::to_string(err) + ": " + msg)
    {
    }

    /** \brief %what() function
     *
     * \returns error message
     */
    const char *what() const noexcept { return m_msg.c_str(); }

private:
    std::string m_msg;
};

/** \brief SQLite Statement handle
 *
 * This object wraps an SQLite statement.  It is iterable in a while
 * loop with the step() function or in a range-based for loop.
 */
class statement
{
    sqlite3 *m_db;
    std::unique_ptr<sqlite3_stmt, std::function<int(sqlite3_stmt *)>> m_stmt;
    std::string m_sql;
    bool m_prepared;
    bool m_valid;

    void handle_error(int err)
    {
        if (err)
            throw error(err, sqlite3_errmsg(m_db));
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
    explicit statement(sqlite3 *db, std::string const &sql)
        : m_db{db}, m_stmt{nullptr, sqlite3_finalize}, m_sql{sql},
          m_prepared{false}, m_valid{true}
    {
    }

    /** \brief Set the query string
     *
     * \param[in] sql   the query string
     * \throws sqlite::error
     */
    void set(std::string const &sql)
    {
        if (m_prepared)
            handle_error(SQLITE_MISUSE);
        m_sql = sql;
    }

    /** \overload void set(std::string const &sql) */
    void set(std::stringstream const &ss) { return set(ss.str()); }

    /** \brief Prepare the statement
     *
     * Prepare the statement for stepping.
     *
     * \throws sqlite::error
     */
    void prepare()
    {
        sqlite3_stmt *pStmt; // is managed below
        auto err = sqlite3_prepare_v2(m_db, m_sql.c_str(), -1, &pStmt, nullptr);
        m_stmt.reset(pStmt);

        handle_error(err);

        m_prepared = true;
    }

    /** \brief Step the statement
     *
     * This will do a busy wait if the statement cannot be stepped
     * immediately.  This is probably a waste of resources but at the
     * same time we do not lose performance in timeouts.
     *
     * \returns true if there is a row, false otherwise
     *
     * \throws sqlite::error
     */
    bool step()
    {
        if (!m_prepared)
            handle_error(SQLITE_MISUSE);
        if (!m_valid)
            handle_error(SQLITE_DONE);

        auto err = SQLITE_BUSY;
        while (err == SQLITE_BUSY)
            err = sqlite3_step(m_stmt.get());

        switch (err) {
        case SQLITE_DONE:
            m_valid = false;
            return false;
        case SQLITE_ROW:
            m_valid = true;
            return true;
        }

        handle_error(err);
        m_valid = false;
        return false;
    }

    /** \brief Reset the statement
     *
     * Reset the statement so you can rebind values.
     *
     * \throws sqlite::error
     */
    void reset()
    {
        handle_error(sqlite3_reset(m_stmt.get()));
        m_valid = true;
        m_prepared = false;
    }

    /** \brief Execute SQLite statements
     *
     * The SQL statements are passed to this function as a string and
     * are executed in-place.
     *
     * \param[in] sql   SQL statements
     * \throws sqlite::error
     */
    void exec(std::string const &sql)
    {
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
    void bind(int where, std::string const &s)
    {
        handle_error(sqlite3_bind_text(m_stmt.get(), where, s.c_str(),
                                       s.length(), SQLITE_STATIC));
    }

    /** \overload void bind(int where, std::string const &s) */
    void bind(int where, std::string &&s)
    {
        handle_error(sqlite3_bind_text(m_stmt.get(), where, s.c_str(),
                                       s.length(), SQLITE_TRANSIENT));
    }

    /** \overload void bind(int where, std::string const &s) */
    void bind(int where, double d)
    {
        handle_error(sqlite3_bind_double(m_stmt.get(), where, d));
    }

    /** \overload void bind(int where, std::string const &s) */
    void bind(int where, int i)
    {
        handle_error(sqlite3_bind_int(m_stmt.get(), where, i));
    }

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
    template <typename T, typename = typename std::enable_if<
                              std::is_floating_point<T>::value>::type>
    double get(int field)
    {
        return sqlite3_column_double(m_stmt.get(), field);
    }

    template <typename T, typename = typename std::enable_if<
                              std::is_same<T, int>::value>::type>
    int get(int field)
    {
        return sqlite3_column_int(m_stmt.get(), field);
    }

    template <typename T, typename = typename std::enable_if<
                              std::is_same<T, std::string>::value>::type>
    std::string get(int field)
    {
        return std::string(reinterpret_cast<char const *>(
            sqlite3_column_text(m_stmt.get(), field)));
    }
#endif

    /** \brief Statement iterator
     *
     * The iterator builds upon Boost Iterator Facade to implement an
     * iterator over the statement with minimal boilerplate.  We refer
     * to the official documentation for the internal functions:
     * http://www.boost.org/libs/iterator/
     */
    class iterator : public boost::iterator_facade<iterator, statement,
                                                   boost::forward_traversal_tag>
    {
    private:
        friend class boost::iterator_core_access;
        statement *m_stmt;
        bool m_done;
        bool m_end;

        void increment()
        {
            if (!m_done)
                m_done = m_stmt->step();
            else
                m_end = true;
        }

        bool equal(iterator const &) const { return m_end; }

        statement &dereference() const
        {
#ifndef NDEBUG
            if (m_end)
                throw std::out_of_range("iterator out of range");
#endif
            return *m_stmt;
        }

    public:
        iterator() : m_stmt{nullptr}, m_done{true}, m_end{true} {}

        explicit iterator(statement *p) : m_stmt{p}, m_done{false}, m_end{false}
        {
            m_stmt->step();
        }
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
class handle final
{
    std::unique_ptr<sqlite3, decltype(&sqlite3_close)> m_db;

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
        : m_db{nullptr, sqlite3_close}
    {
        sqlite3 *tmp_db;
        auto err = sqlite3_open_v2(filename.c_str(), &tmp_db, flags, nullptr);
        m_db.reset(tmp_db);

        if (err)
            throw error(err, sqlite3_errmsg(m_db.get()));
    }

    /** \brief Execute SQLite statements
     *
     * The SQL statements are passed to this function as a string and
     * are executed in-place.
     *
     * \deprecated This is deprecated and should not be used anymore.
     * The usage of prepared statements in encouraged instead.  Calls
     * to this function are equivalent to and can hence be replaced by
     * \code
     * sqlite::statement stmt(db);
     * stmt.exec(sql);
     * \endcode
     *
     * \param[in] sql   SQL statements
     */
#if defined(_MSC_VER)
    __declspec(deprecated)
#else
    __attribute__((deprecated))
#endif
        void exec(std::string const &sql)
    {
        auto err = sqlite3_exec(*this, sql.c_str(), nullptr, nullptr, nullptr);
        if (err)
            throw error(err, sqlite3_errmsg(m_db.get()));
    }
};

} // namespace sqlite

#endif // SQLITE_BINDINGS
