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

class SQLite3Result {
public:
    class iterator {
    public:
        iterator(const SQLite3Result* result, int pos, int nColumn);
        bool operator!= (const iterator& other) const;
        std::unique_ptr<std::stringstream> operator* () const;
        const iterator& operator++ ();
    private:
        const SQLite3Result* result;
        int pos, nColumn;
        int rc;
    };
    SQLite3Result();
    ~SQLite3Result();
    size_t size() const;
    SQLite3Result::iterator begin() const;
    SQLite3Result::iterator end() const;
    std::unique_ptr<std::stringstream> first() const;
    std::unique_ptr<std::stringstream> header() const;
    std::unique_ptr<std::stringstream> getRow(int pos) const;
    char **azResult;
    int nRow, nColumn;
};

class SQLite3 {
public:
    SQLite3(const std::string filename);
    SQLite3(const std::string filename, int flags);
    ~SQLite3();
    SQLite3Result query(const std::string sql);
    void exec(const std::string sql);
private:
    sqlite3 *db;
    char *zErrMsg;
    int rc;
};

#endif // SQLITE_BINDINGS
