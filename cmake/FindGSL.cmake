# Find the GNU Scientific Library
#
# GSL_LIBRARIES      - gsl and gslcblas libraries
# GSL_INCLUDE_DIRS   - directory with all header files

# Include these modules to handle the QUIETLY and REQUIRED arguments.
include (FindPackageHandleStandardArgs)

find_package(PkgConfig)
pkg_check_modules(GSL gsl)

find_path( GSL_INCLUDE_DIR
  NAMES gsl/gsl_sf.h
  HINTS ${GSL_INCLUDEDIR} )

find_library( GSL_LIBRARY
  NAMES gsl
  HINTS ${GSL_LIBDIR} )

find_library( GSL_CBLAS_LIBRARY
  NAMES gslcblas cblas
  HINTS ${GSL_LIBDIR} )

set( GSL_INCLUDE_DIRS ${GSL_INCLUDE_DIR} )
set( GSL_LIBRARIES ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY} )

# Set the FOUND variable to TRUE if all listed variables are set.
find_package_handle_standard_args( GSL DEFAULT_MSG GSL_LIBRARIES GSL_INCLUDE_DIRS )

mark_as_advanced( GSL_LIBRARIES GSL_INCLUDE_DIRS )
