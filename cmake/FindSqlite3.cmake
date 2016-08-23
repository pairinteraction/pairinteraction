# Find SQLite3
#
# SQLITE3_LIBRARY       - the sqlite3 library
# SQLITE3_INCLUDE_DIR   - path including sqlite3.h
# SQLITE3_BINARY        - sqlite3 executable

# Include these modules to handle the QUIETLY and REQUIRED arguments.
include (FindPackageHandleStandardArgs)

find_path( SQLITE3_INCLUDE_DIR NAMES sqlite3.h )

find_library( SQLITE3_LIBRARY NAMES sqlite3 )

find_program( SQLITE3_BINARY NAMES sqlite3 )

# Set the FOUND variable to TRUE if all listed variables are set.
find_package_handle_standard_args( SQLITE3 DEFAULT_MSG SQLITE3_LIBRARY SQLITE3_INCLUDE_DIR SQLITE3_BINARY )

mark_as_advanced( SQLITE3_LIBRARY SQLITE3_INCLUDE_DIR SQLITE3_BINARY )
