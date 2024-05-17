find_package(PkgConfig QUIET)
pkg_search_module(httplib QUIET IMPORTED_TARGET cpp-httplib)
if ("-DCPPHTTPLIB_OPENSSL_SUPPORT" IN_LIST httplib_CFLAGS_OTHER)
    set(httplib_OpenSSL_FOUND TRUE)
endif()
if ("-DCPPHTTPLIB_ZLIB_SUPPORT" IN_LIST httplib_CFLAGS_OTHER)
    set(httplib_ZLIB_FOUND TRUE)
endif()
if ("-DCPPHTTPLIB_BROTLI_SUPPORT" IN_LIST httplib_CFLAGS_OTHER)
    set(httplib_Brotli_FOUND TRUE)
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(httplib
    REQUIRED_VARS httplib_FOUND
    VERSION_VAR httplib_VERSION
    HANDLE_COMPONENTS
)

add_library(httplib::httplib ALIAS PkgConfig::httplib)
