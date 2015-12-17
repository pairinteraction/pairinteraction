#ifndef SQLITE_BINDINGS
#define SQLITE_BINDINGS

#include <iostream>
#include <sstream>
#include <sqlite3.h>
#include "Iter.h"

class SQLite3Result {
public:
  char **azResult;
  int nRow, nColumn;
  SQLite3Result() : azResult(NULL), nRow(0), nColumn(0) {}
  ~SQLite3Result() { if (azResult) sqlite3_free_table(azResult); }
};

std::ostream& operator<<(std::ostream& os, const SQLite3Result& res);


class SQLite3Row {
public:
    SQLite3Row(char **azResult, int pos, int nColumn) : azResult(azResult), pos(pos), nColumn(nColumn) {
    }
    std::stringstream operator[](int index) const {
        std::stringstream ss;
        ss << azResult[(pos+1)*nColumn+index];
        return ss;
    }
private:
    char **azResult;
    int pos, nColumn;
};

class SQLite3Iter {
public:
    SQLite3Iter(char **azResult, int pos, int nColumn) : azResult(azResult), pos(pos), nColumn(nColumn) {
    }
    bool operator!= (const SQLite3Iter& other) const {
        return pos != other.pos;
    }
    const SQLite3Row operator* () const {
        return SQLite3Row(azResult, pos, nColumn);
    }
    const SQLite3Iter& operator++ () {
        ++pos;
        return *this;
    }
private:
    char **azResult;
    int pos, nColumn;
    int rc;
};

class SQLite3Result2 {
public:
    SQLite3Result2() : azResult(NULL), nRow(0), nColumn(0) {
    }
    ~SQLite3Result2() {
        sqlite3_free_table(azResult);
    }
    size_t size() const {
        return nRow;
    }
    SQLite3Iter begin() const {
        return SQLite3Iter(azResult, 0, nColumn);
    }
    SQLite3Iter end() const {
        return SQLite3Iter(azResult, nRow, nColumn);
    }
    SQLite3Row first() {
        return SQLite3Row(azResult, 0, nColumn);
    }
    SQLite3Row header() {
        return SQLite3Row(azResult, -1, nColumn);
    }
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
  SQLite3Result2 query2(const char* sql);
  void exec(const char* sql);
};

#endif // SQLITE_BINDINGS
