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
     * \param[in] msg   error message
     */
    explicit error(std::string const &msg)
        : m_msg(std::string("SQLite error: ") + msg)
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
 * loop with the step() function.
 */
class statement
{
    sqlite3 *m_db;
    std::unique_ptr<sqlite3_stmt, std::function<int(sqlite3_stmt *)>> m_stmt;
    std::string m_sql;
    bool m_prepared;

    void handle_error(int err)
    {
        if (err)
            throw error(sqlite3_errmsg(m_db));
    }

public:
    /** \brief Constructor
     *
     * This is called in the query function of the the handle object.
     *
     * \param[in] db    unmanaged raw pointer to the sqlite database
     */
    explicit statement(sqlite3 *db)
        : m_db{db}, m_stmt{nullptr, sqlite3_finalize}, m_sql{}, m_prepared{
                                                                    false}
    {
    }

    /** \brief Constructor
     *
     * This is called in the query function of the the handle object.
     * This variant also sets a default query.
     *
     * \param[in] db    unmanaged raw pointer to the sqlite database
     * \param[in] sql   an initial query as a string
     */
    explicit statement(sqlite3 *db, std::string const &sql)
        : m_db{db}, m_stmt{nullptr, sqlite3_finalize}, m_sql{sql}, m_prepared{
                                                                       false}
    {
    }

    /** \brief Set the query string
     *
     * \param[in] sql   the query string
     */
    void set(std::string const &sql)
    {
        if (m_prepared)
            throw error("Cannot set query on a prepared statement!");
        m_sql = sql;
    }

    /** \overload void set(std::string const &sql) */
    void set(std::stringstream const &ss) { return set(ss.str()); }

    /** \brief Prepare the statement
     *
     * Prepare the statement for stepping.
     */
    void prepare()
    {
        m_prepared = true;

        char const *zTail;   // only points into m_sql, no management needed
        sqlite3_stmt *pStmt; // is managed below

        auto err = sqlite3_prepare_v2(m_db, m_sql.c_str(), m_sql.length(),
                                      &pStmt, &zTail);
        m_stmt.reset(pStmt);

        handle_error(err);
    }

    /** \brief Step the statement
     *
     * This will do a busy wait if the statement cannot be stepped
     * immediately.  This is probably a waste of resources but at the
     * same time we do not lose performance in timeouts.
     */
    bool step()
    {
        if (!m_prepared)
            throw error("Cannot step an unprepared statement!");

        auto err = SQLITE_BUSY;
        while (err == SQLITE_BUSY)
            err = sqlite3_step(m_stmt.get());

        switch (err) {
        case SQLITE_DONE:
            return false;
        case SQLITE_ROW:
            return true;
        }

        handle_error(err);
        return false;
    }

    /** \brief Reset the statement
     *
     * Reset the statement so you can rebind values.
     */
    void reset()
    {
        handle_error(sqlite3_reset(m_stmt.get()));
        m_prepared = false;
    }

    /** \brief Bind a value to a position in the query
     *
     * \param[in] where   position in the query (one-based)
     * \param[in] s       string to bind
     */
    void bind(int where, std::string const &s)
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
};

/** \brief SQLite Database handle
 *
 * This object handles the connection to the underlying SQLite databse.  It
 * provides a high-level object oriented RAII interface to the C-bindings of
 * SQLite3.
 */
class handle final
{
    typedef std::unique_ptr<sqlite3, std::function<int(sqlite3 *)>> sqlite3_ptr;

public:
    /** \brief Conversion operator
     *
     * The handle object implicitly behaves like the underlying
     * pointer.
     */
    operator sqlite3 *() { return db.get(); }

    /** \brief Constructor
     *
     * This is the canonical constructor.  It takes a filename and opens a
     * database connection with the default parameters.
     *
     * \param[in] filename    fully-qualified filename of the database
     */
    explicit handle(std::string const &filename) : db(nullptr, sqlite3_close)
    {
        sqlite3 *_db;
        auto err = sqlite3_open(filename.c_str(), &_db);
        db = sqlite3_ptr(_db, sqlite3_close);

        if (err)
            throw error(sqlite3_errmsg(db.get()));
    }

    /** \brief Extended Constructor
     *
     * This is the alternative constructor.  It takes a filename and a bitmask
     * of parameters to open the database connection.  This is useful if you
     * have to open the connection, e.g., read-only.
     *
     * \param[in] filename    fully-qualified filename of the database
     * \param[in] flags       options to open the database (bitmask)
     * \throws sqlite::error
     */
    explicit handle(std::string const &filename, int flags)
        : db(nullptr, sqlite3_close)
    {
        sqlite3 *_db;
        auto err = sqlite3_open_v2(filename.c_str(), &_db, flags, nullptr);
        db = sqlite3_ptr(_db, sqlite3_close);

        if (err) {
            throw error(sqlite3_errmsg(db.get()));
        }
    }

    /** \brief Execute SQLite statements
     *
     * The SQL statements are passed to this function as a string (or
     * stringstream) and are executed in-place.
     *
     * \param[in] sql   SQL statements
     */
    void exec(std::string const &sql)
    {
        statement stmt(db.get(), sql);
        stmt.prepare();
        stmt.step();
    }

    /** \overload void exec(std::string const &sql) */
    void exec(std::stringstream const &ss) { return exec(ss.str()); }

private:
    sqlite3_ptr db;
};

} // namespace sqlite

#endif // SQLITE_BINDINGS
