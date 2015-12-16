#include <iostream>
#include <sstream>
#include "SQLite.hpp"

std::ostream& operator<<(std::ostream& os, const SQLite3Result& res) {
  os << "Result has " << res.nRow << " rows and "
     << res.nColumn << " columns" << std::endl;
  for (int i = 0; i <= res.nRow*res.nColumn; i += res.nColumn) {
    for (int j = 0; j < res.nColumn; j++)
      os << res.azResult[i+j] << (j < res.nColumn-1 ? '|' : '\n');
  }
  return os;
}

SQLite3Result SQLite3::query(const char* sql) {
  SQLite3Result res;
  rc = sqlite3_get_table(db, sql, &res.azResult, &res.nRow, &res.nColumn, &zErrMsg);
  if ( rc != SQLITE_OK ) {
    std::cerr << "SQL error: " << zErrMsg << std::endl;
    sqlite3_free(zErrMsg);
  }
  return res;
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

SQLite3Result2 SQLite3::query2(const char* sql) {
  SQLite3Result2 result;
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
