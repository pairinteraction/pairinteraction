#include <iostream>
#include <sstream>
#include "SQLite.hpp"

SQLite3Result::iterator::iterator(const SQLite3Result* result, int pos, int nColumn) : result(result), pos(pos), nColumn(nColumn) {
}
bool SQLite3Result::iterator::operator!= (const SQLite3Result::iterator& other) const {
    return pos != other.pos;
}
std::unique_ptr<std::stringstream> SQLite3Result::iterator::operator* () const {
    return result->getRow(pos);
}
const SQLite3Result::iterator& SQLite3Result::iterator::operator++ () {
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
SQLite3Result::iterator SQLite3Result::begin() const {
    return SQLite3Result::iterator(this, 0, nColumn);
}
SQLite3Result::iterator SQLite3Result::end() const {
    return SQLite3Result::iterator(this, nRow, nColumn);
}
std::unique_ptr<std::stringstream> SQLite3Result::first() const {
    return getRow(0);
}
std::unique_ptr<std::stringstream> SQLite3Result::header() const {
    return getRow(-1);
}
std::unique_ptr<std::stringstream> SQLite3Result::getRow(int pos) const {
    std::stringstream ss;
    std::string spacer = "";
    for (int i = 0; i < nColumn; ++i) {
        ss << spacer << azResult[(pos+1)*nColumn+i];
        spacer = " ";
    }
    return std::unique_ptr<std::stringstream>(new std::stringstream(ss.str()));
}

SQLite3::SQLite3(const std::string filename) : zErrMsg(NULL) {
    rc = sqlite3_open(filename.c_str(), &db);
    if ( rc ) {
        std::cerr << "Can't open database: " << sqlite3_errmsg(db) << std::endl;
    }
}
SQLite3::~SQLite3() {
    sqlite3_close(db);
}
SQLite3Result SQLite3::query(const std::string sql) {
    SQLite3Result result;
    rc = sqlite3_get_table(db, sql.c_str(), &result.azResult, &result.nRow, &result.nColumn, &zErrMsg);
    if ( rc != SQLITE_OK ) {
        std::cerr << "SQL error: " << zErrMsg << std::endl;
        sqlite3_free(zErrMsg);
    }
    return result;
}
void SQLite3::exec(const std::string sql) {
    rc = sqlite3_exec(db, sql.c_str(), NULL, NULL, &zErrMsg);
    if ( rc != SQLITE_OK ) {
        std::cerr << "SQL error: " << zErrMsg << std::endl;
        sqlite3_free(zErrMsg);
    }
}
