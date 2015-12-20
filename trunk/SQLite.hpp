#ifndef SQLITE_BINDINGS
#define SQLITE_BINDINGS

#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <sqlite3.h>
#include "Iter.h"

class SQLite3Result {
public:
    class iterator {
    public:
        iterator(const SQLite3Result* result, int pos, int nColumn);
        bool operator!= (const iterator& other) const;
        std::stringstream operator* () const;
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
    std::stringstream first() const;
    std::stringstream header() const;
    std::stringstream getRow(int pos) const;
    char **azResult;
    int nRow, nColumn;
};

class SQLite3 {
public:
    SQLite3(const std::string filename);
    ~SQLite3();
    SQLite3Result query(const std::string sql);
    void exec(const std::string sql);
private:
    sqlite3 *db;
    char *zErrMsg;
    int rc;
};

#endif // SQLITE_BINDINGS
