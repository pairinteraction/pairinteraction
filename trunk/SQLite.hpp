#ifndef SQLITE_BINDINGS
#define SQLITE_BINDINGS

#include <iostream>
#include <sstream>
#include <vector>
#include <sqlite3.h>
#include "Iter.h"

class SQLite3Result;

class SQLite3Iter {
public:
    SQLite3Iter(const SQLite3Result* result, int pos, int nColumn);
    bool operator!= (const SQLite3Iter& other) const;
    std::stringstream operator* () const;
    const SQLite3Iter& operator++ ();
private:
    const SQLite3Result* result;
    int pos, nColumn;
    int rc;
};

class SQLite3Result {
public:
    SQLite3Result();
    ~SQLite3Result();
    size_t size() const;
    SQLite3Iter begin() const;
    SQLite3Iter end() const;
    std::stringstream first() const;
    std::stringstream header() const;
    std::stringstream getRow(int pos) const;
    char **azResult;
    int nRow, nColumn;
};

class SQLite3 {
private:
    sqlite3 *db;
    char *zErrMsg;
    int rc;
public:
    SQLite3(const char* filename);
    ~SQLite3();
    SQLite3Result query(const char* sql);
    void exec(const char* sql);
};

#endif // SQLITE_BINDINGS
