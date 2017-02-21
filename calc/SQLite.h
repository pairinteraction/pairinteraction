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

#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <sqlite3.h>
#include <memory>
#include "Iter.h"


namespace sqlite {

/** \brief SQLite error
 *
 * Exception to identify SQLite errors.
 * \code
 * try
 * {
 *   // code with sqlite
 * }
 * catch (sqlite::sqlite_error &e)
 * {
 *   // handle exception
 * }
 * \endcode
 */
struct sqlite_error : public std::exception {
public:
    /** \brief Constructor
     *
     * \param[in] msg   error message
     */
    explicit sqlite_error(std::string const& msg)
        : m_msg(std::string("SQLite error: ") + msg)
    {}

    /** \brief %what() function
     *
     * \returns error message
     */
    const char* what() const noexcept
    {
        return m_msg.c_str();
    }

private:
    std::string m_msg;
};


/** \brief Result table
 *
 * This object wraps the result table.  It is iterable and can be used in
 * range-based for loops.
 */
class result {
    friend class handle;
public:
    /** \brief Row of the table
     *
     * The rows of the table are saved as a string inside a stringstream and
     * deserialized using the >> operator.
     */
    class row {
    public:
        /** \brief Constructor
         *
         * Intialize a row from a string.
         *
         * \param[in] s   contents of the row
         */
        explicit row(std::string const& s)
            : m_row( new std::stringstream(s) )
        {}

        /** \brief Conversion operator
         *
         * A row can behave like a string.
         */
        operator std::string() const
        {
            return m_row->str();
        }

        /** \brief Deserialization operator
         *
         * Single elements can be piped out.
         *
         * \param[out] other    target to write the element
         * \returns the row itself
         */
        template < typename T >
        row& operator>> ( T& other )
        {
            *m_row >> other;
            return *this;
        }

        /** \brief Deserialization operator */
        template < typename T >
        row const& operator>> ( T& other ) const
        {
            *m_row >> other;
            return *this;
        }

    private:
        std::unique_ptr < std::stringstream > m_row;
    }; // class row

    /** \brief Table iterator
     *
     * The iterator gives access to the rows of the table.
     */
    class iterator {
    public:
        /** \brief Constructor
         *
         * Construct an iterator for a result table.
         *
         * \param[in] res       pointer to the result table
         * \param[in] pos       position of the row
         * \param[in] nColumn   number of columns
         */
        explicit iterator(result const* res, int pos, int nColumn)
            : res(res), pos(pos), nColumn(nColumn)
        {}

        /** \brief Comparison operator
         *
         * Required to find the end of the range in a range-based loop.
         *
         * \param[in] other   iterator to compare to
         */
        bool operator!= (const iterator& other) const
        {
            return pos != other.pos;
        }

        /** \brief Dereference operator
         *
         * Get the current row.
         *
         * \returns The current row.
         */
        row operator* () const
        {
            return res->getRow(pos);
        }

        /** \brief Increment operator
         *
         * Advance iterator to the next row.
         *
         * \returns The advanced iterator
         */
        iterator& operator++ ()
        {
            ++pos;
            return *this;
        }

    private:
        result const* res;
        int pos, nColumn;
    }; // class iterator

    /** \brief Constructor
     *
     * Construct the result table of a query.
     */
    explicit result( char **azResult, int nRow, int nColumn )
        : azResult(azResult, sqlite3_free_table)
        , nRow(nRow)
        , nColumn(nColumn)
    {}

    /** \brief Number of rows
     *
     * \returns Number of rows
     */
    size_t size() const
    {
        return nRow;
    }

    /** \brief Iterator pointing to the first row
     *
     * \returns Iterator pointing to the first row
     */
    result::iterator begin() const
    {
        return result::iterator(this, 0, nColumn);
    }

    /** \brief Iterator pointing to the last row
     *
     * \returns Iterator pointing to the last row
     */
    result::iterator end() const
    {
        return result::iterator(this, nRow, nColumn);
    }

    /** \brief First row
     *
     * \returns First row
     */
    result::row first() const
    {
        return getRow(0);
    }

    /** \brief Header row
     *
     * \returns Header row
     */
    result::row header() const
    {
        return getRow(-1);
    }

    /** \brief Get a specific row
     *
     * \param[in] pos   position of the row
     * \returns Specified row
     * \throws sqlite_error
     */
    result::row getRow(int pos) const
    {
        int p = pos + 1;
        if ( p > nRow || p < 0 )
        {
            throw sqlite_error("Position index out of range");
        }
        std::string output;
        std::string spacer = "";
        for (int i = 0; i < nColumn; ++i) {
            output += spacer + azResult.get()[p*nColumn+i];
            spacer = " ";
        }
        return result::row(output);
    }

private:
    std::unique_ptr < char*, std::function< void(char**) > > azResult;
    int nRow, nColumn;
};


/** \brief SQLite Database handle
 *
 * This object handles the connection to the underlying SQLite databse.  It
 * provides a high-level object oriented RAII interface to the C-bindings of
 * SQLite3.
 */
class handle {
    typedef std::unique_ptr < sqlite3, std::function<int(sqlite3*)> > sqlite3_ptr;
public:
    /** \brief Constructor
     *
     * This is the canonical constructor.  It takes a filename and opens a
     * database connection with the default parameters.
     *
     * \param[in] filename    fully-qualified filename of the database
     */
    explicit handle(std::string const& filename)
        : db(nullptr, sqlite3_close)
    {
        sqlite3* _db;
        auto err = sqlite3_open(filename.c_str(), &_db);
        db = sqlite3_ptr(_db, sqlite3_close);

        if ( err )
        {
            throw sqlite_error( sqlite3_errmsg( db.get() ) );
        }
    }

    /** \brief Extended Constructor
     *
     * This is the alternative constructor.  It takes a filename and a bitmask
     * of parameters to open the database connection.  This is useful if you
     * have to open the connection, e.g., read-only.
     *
     * \param[in] filename    fully-qualified filename of the database
     * \param[in] flags       options to open the database (bitmask)
     * \throws sqlite_error
     */
    explicit handle(std::string const& filename, int flags)
        : db(nullptr, sqlite3_close)
    {
        sqlite3* _db;
        auto err = sqlite3_open_v2(filename.c_str(), &_db, flags, nullptr);
        db = sqlite3_ptr(_db, sqlite3_close);

        if ( err )
        {
            throw sqlite_error( sqlite3_errmsg( db.get() ) );
        }
    }

    /** \brief Place an SQLite query
     *
     * The SQL statements are passed to this function as a string (or
     * stringstream) and are executed in-place.  The string has to be prepared
     * before passing to this function.  The table containing the result is
     * returned in an sqlite::result object.
     *
     * \param[in] sql   SQL statements
     * \returns result of the query
     * \throws sqlite_error
     */
    result query(std::string const& sql)
    {
        char *zErrMsg;
        char **azResult;
        int nRow, nColumn;
        auto err = sqlite3_get_table( db.get(),    // open database handle
                                      sql.c_str(), // sql query
                                      &azResult,   // result table
                                      &nRow,       // number of rows
                                      &nColumn,    // number of columns
                                      &zErrMsg );  // error message

        result res( azResult, nRow, nColumn );

        if ( err != SQLITE_OK )
        {
            std::string msg(zErrMsg);
            sqlite3_free(zErrMsg);
            throw sqlite_error(msg);
        }

        return res;
    }

    /** \brief Place an SQLite query */
    result query(std::stringstream const& ss)
    {
        return query(ss.str());
    }

    /** \brief Execute SQLite statements
     *
     * The SQL statements are passed to this function as a string (or
     * stringstream) and are executed in-place.  The string has to be prepared
     * before passing to this function.  This function discards the result table.
     *
     * \param[in] sql   SQL statements
     * \throws sqlite_error
     */
    void exec(std::string const& sql)
    {
        char *zErrMsg;
        auto err = sqlite3_exec( db.get(),    // open database handle
                                 sql.c_str(), // sql query
                                 nullptr,     // callback
                                 nullptr,     // 1st argument to callback
                                 &zErrMsg );  // error message

        if ( err != SQLITE_OK )
        {
            std::string msg(zErrMsg);
            sqlite3_free(zErrMsg);
            throw sqlite_error(msg);
        }
    }

    /** \brief Execute SQLite statements */
    void exec(std::stringstream const& ss)
    {
        return exec(ss.str());
    }

private:
    sqlite3_ptr db;
};

} // namespace sqlite


#endif // SQLITE_BINDINGS
