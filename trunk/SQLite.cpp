#include <iostream>
#include <sstream>
#include "SQLite.hpp"

SQLite3Iter::SQLite3Iter(const SQLite3Result* result, int pos, int nColumn) : result(result), pos(pos), nColumn(nColumn) {
}
bool SQLite3Iter::operator!= (const SQLite3Iter& other) const {
    return pos != other.pos;
}
std::stringstream SQLite3Iter::operator* () const {
    return result->getRow(pos);
}
const SQLite3Iter& SQLite3Iter::operator++ () {
    ++pos;
    return *this;
}

SQLite3Result::SQLite3Result() : azResult(NULL), nRow(0), nColumn(0) {
}
SQLite3Result::~SQLite3Result() {
    sqlite3_free_table(azResult);
}
size_t SQLite3Result::size() const {
    return nRow;
}
SQLite3Iter SQLite3Result::begin() const {
    return SQLite3Iter(this, 0, nColumn);
}
SQLite3Iter SQLite3Result::end() const {
    return SQLite3Iter(this, nRow, nColumn);
}
std::stringstream SQLite3Result::first() const {
    return getRow(0);
}
std::stringstream SQLite3Result::header() const {
    return getRow(1);
}
std::stringstream SQLite3Result::getRow(int pos) const {
    std::stringstream ss;
    for (int i = 0; i < nColumn; ++i) {
        ss << azResult[(pos+1)*nColumn+i]<< " ";
    }
    return ss;
}

SQLite3::SQLite3(const char* filename) : zErrMsg(NULL) {
    rc = sqlite3_open(filename, &db);
    if ( rc ) {
        std::cerr << "Can't open database: " << sqlite3_errmsg(db) << std::endl;
    }
}
SQLite3::~SQLite3() {
    sqlite3_close(db);
}
SQLite3Result SQLite3::query(const char* sql) {
    SQLite3Result result;
    rc = sqlite3_get_table(db, sql, &result.azResult, &result.nRow, &result.nColumn, &zErrMsg);
    if ( rc != SQLITE_OK ) {
        std::cerr << "SQL error: " << zErrMsg << std::endl;
        sqlite3_free(zErrMsg);
    }
    return result;
}
void SQLite3::exec(const char* sql) {
    rc = sqlite3_exec(db, sql, NULL, NULL, &zErrMsg);
    if ( rc != SQLITE_OK ) {
        std::cerr << "SQL error: " << zErrMsg << std::endl;
        sqlite3_free(zErrMsg);
    }
}
