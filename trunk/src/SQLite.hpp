#ifndef SQLITE_BINDINGS
#define SQLITE_BINDINGS

#include <iostream>
#include <sstream>
#include <sqlite3.h>

class SQLite3Result {
public:
  char **azResult;
  int nRow, nColumn;
  SQLite3Result() : azResult(NULL), nRow(0), nColumn(0) {};
  ~SQLite3Result() { if (azResult) sqlite3_free_table(azResult); };
};


std::ostream& operator<<(std::ostream& os, const SQLite3Result& res);


class SQLite3 {
private:
  sqlite3 *db;
  char *zErrMsg;
  int rc;
public:
  SQLite3(const char* filename);
  ~SQLite3();
  SQLite3Result query(const char* sql);
};

#endif // SQLITE_BINDINGS
