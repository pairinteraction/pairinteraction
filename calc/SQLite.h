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

struct sqlite_error : public std::exception {
public:
    explicit sqlite_error(std::string msg)
        : m_msg(std::string("SQLite error: ") + msg)
    {}

    const char* what() const noexcept
    {
        return m_msg.c_str();
    }

private:
    std::string m_msg;
};


class result {
public:
    class iterator {
    public:
        explicit iterator(const result* res, int pos, int nColumn)
            : res(res), pos(pos), nColumn(nColumn)
        {}

        bool operator!= (const iterator& other) const
        {
            return pos != other.pos;
        }

        std::unique_ptr<std::stringstream> operator* () const
        {
            return res->getRow(pos);
        }

        const iterator& operator++ ()
        {
            ++pos;
            return *this;
        }

    private:
        const result* res;
        int pos, nColumn;
        int rc;
    };

    explicit result()
        : azResult(NULL), nRow(0), nColumn(0)
    {}

    result( result& r )
    {
        swap( r );
        // invalidate copied result to prevent deletion
        r.azResult = nullptr;
    }

    result( result&& r )
    {
        swap( r );
        // invalidate copied result to prevent deletion
        r.azResult = nullptr;
    }

    void swap( result& r )
    {
        std::swap ( nRow    , r.nRow     );
        std::swap ( nColumn , r.nColumn  );
        std::swap ( azResult, r.azResult );
    }

    void swap( result&& r )
    {
        std::swap ( nRow    , r.nRow     );
        std::swap ( nColumn , r.nColumn  );
        std::swap ( azResult, r.azResult );
    }

    ~result()
    {
        sqlite3_free_table(azResult);
    }

    size_t size() const
    {
        return nRow;
    }

    result::iterator begin() const
    {
        return result::iterator(this, 0, nColumn);
    }

    result::iterator end() const
    {
        return result::iterator(this, nRow, nColumn);
    }

    std::unique_ptr<std::stringstream> first() const
    {
        return getRow(0);
    }

    std::unique_ptr<std::stringstream> header() const
    {
        return getRow(-1);
    }

    std::unique_ptr<std::stringstream> getRow(int pos) const
    {
        std::unique_ptr<std::stringstream> ss( new std::stringstream() );
        std::string spacer = "";
        for (int i = 0; i < nColumn; ++i) {
            *ss << spacer << azResult[(pos+1)*nColumn+i];
            spacer = " ";
        }
        return ss;
    }

    char **azResult;
    int nRow, nColumn;

private:
    result& operator=( result const& );
    result& operator=( result&& );
};

class handle {
public:
    explicit handle(const std::string filename)
        : zErrMsg(NULL)
    {
        if ( sqlite3_open(filename.c_str(), &db) )
        {
            sqlite3_close(db);
            throw sqlite_error( sqlite3_errmsg(db) );
        }
    }

    explicit handle(const std::string filename, int flags)
        : zErrMsg(NULL)
    {
        if ( sqlite3_open_v2(filename.c_str(), &db, flags, NULL) )
        {
            sqlite3_close(db);
            throw sqlite_error( sqlite3_errmsg(db) );
        }
    }

    ~handle()
    {
        sqlite3_close(db);
    }

    result query(const std::string sql)
    {
        result res;
        if ( sqlite3_get_table(db, sql.c_str(), &res.azResult, &res.nRow, &res.nColumn, &zErrMsg) != SQLITE_OK ) {
            std::string msg("SQL error: ");
            msg += zErrMsg;
            sqlite3_free(zErrMsg);
            throw sqlite_error(msg);
        }
        return res;
    }

    void exec(const std::string sql)
    {
        if ( sqlite3_exec(db, sql.c_str(), NULL, NULL, &zErrMsg) != SQLITE_OK ) {
            std::string msg("SQL error: ");
            msg += zErrMsg;
            sqlite3_free(zErrMsg);
            throw sqlite_error(msg);
        }
    }

private:
    handle( handle const& );
    handle& operator=( handle const& );

    handle( handle&& );
    handle& operator=( handle&& );

    sqlite3 *db;
    char *zErrMsg;
    int rc;
};

} // namespace sqlite


#endif // SQLITE_BINDINGS
