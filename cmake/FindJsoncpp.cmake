# Find JsonCPP
#
# JSONCPP_LIBRARIES      - the jsoncpp library
# JSONCPP_INCLUDE_DIRS   - directory with all header files

# Include these modules to handle the QUIETLY and REQUIRED arguments.
include (FindPackageHandleStandardArgs)

find_path( JSONCPP_INCLUDE_DIRS NAMES json/json.h PATH_SUFFIXES jsoncpp )

find_library( JSONCPP_LIBRARIES NAMES jsoncpp )

# Set the FOUND variable to TRUE if all listed variables are set.
find_package_handle_standard_args( JSONCPP DEFAULT_MSG JSONCPP_LIBRARIES JSONCPP_INCLUDE_DIRS )

mark_as_advanced( JSONCPP_LIBRARIES JSONCPP_INCLUDE_DIRS )
