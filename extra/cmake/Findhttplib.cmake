include(FindPackageHandleStandardArgs)

find_package(httplib QUIET CONFIG)
if (httplib_CONSIDERED_CONFIGS)
    find_package_handle_standard_args(httplib
        HANDLE_COMPONENTS
        CONFIG_MODE
    )
    if (httplib_VERSION VERSION_GREATER_EQUAL 0.11.0)
        set_property(TARGET httplib::httplib APPEND PROPERTY INTERFACE_COMPILE_DEFINITIONS HTTPLIB_USES_STD_STRING)
    endif()
else()
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
    find_package_handle_standard_args(httplib
        VERSION_VAR httplib_VERSION
        HANDLE_COMPONENTS
    )
    if (httplib_VERSION VERSION_GREATER_EQUAL 0.11.0)
        set_property(TARGET PkgConfig::httplib APPEND PROPERTY INTERFACE_COMPILE_DEFINITIONS HTTPLIB_USES_STD_STRING)
    endif()
    if (httplib_FOUND)
        add_library(httplib::httplib ALIAS PkgConfig::httplib)
    endif()
endif()
