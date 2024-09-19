%global forgeurl0 https://github.com/yhirose/cpp-httplib
%undefine __cmake_in_source_build

%bcond_without tests
%bcond_with    online
# Compiled version in shared library.
# Does not have any so-version, therefore not default
%bcond_with compile

%if %{without compile}
%undefine __cmake_in_source_build
%global debug_package %{nil}
%endif

Name:           cpp-httplib
Version:        0.12.5
Release:        1%{?dist}

Summary:        A C++11 single-file header-only cross platform HTTP/HTTPS library
License:        MIT
URL:            https://github.com/yhirose/cpp-httplib
VCS:            git:%{forgeurl0}
Source0:        %{url}/archive/v%{version}/%{name}-%{version}.tar.gz

BuildRequires:  redhat-rpm-config
BuildRequires:  gcc-c++
BuildRequires:  cmake3
BuildRequires:  openssl-devel
BuildRequires:  zlib-devel
BuildRequires:  brotli-devel
%if %{with tests}
BuildRequires:  openssl
BuildRequires:  gtest-devel
%endif

%description
A C++11 single-file header-only cross platform HTTP/HTTPS library.

It's extremely easy to setup. Just include the httplib.h file in your code!

%package devel
Summary:        A C++11 single-file header-only cross platform HTTP/HTTPS library
Requires:       openssl-devel%{?_isa} zlib-devel%{?_isa} brotli-devel%{?_isa}
%if %{with compile}
Requires:       %{name}%{?_isa} = %{version}-%{release}
%else
Provides:       %{name}-static = %{version}-%{release}
%endif

%description devel
A C++11 single-file header-only cross platform HTTP/HTTPS library.

It's extremely easy to setup. Just include the httplib.h file in your code!

NOTE: This is a multi-threaded 'blocking' HTTP library.
If you are looking for a 'non-blocking' library, this is not the one that you want.

Development files only.

%prep
%autosetup -p1


%build
%cmake3 \
%if %{with compile}
    -DBUILD_SHARED_LIBS=ON -DHTTPLIB_COMPILE=ON \
%endif
%if %{with tests}
    -DHTTPLIB_TEST=ON \
%endif
#
%cmake3_build


%install
%cmake3_install


%check
%if %{with tests}
# multiple threads fails many tests
EXCLUDE=""
%if %{without online}
  EXCLUDE+=" --exclude-regex _Online$"
%endif
%if 0%{?rhel} && 0%{?rhel} < 9
%ifarch s390x x86_64
  # fails reliably on this architecture
  EXCLUDE+="|ServerTest.Gzip|ServerTest.PutLargeFileWithGzip2$"
%endif
%endif
%ctest3 --parallel 1 $EXCLUDE
%endif


%if %{with compile}
%files
%license LICENSE
%doc README.md
# TODO: should use so-versioned library here, but upstream
# prefers header-only mode.
%{_libdir}/libhttplib.so.*
%endif

%files devel
%if %{without compile}
%license LICENSE
%doc README.md
%else
%{_libdir}/libhttplib.so
%endif
%{_includedir}/httplib.h
%{_libdir}/cmake/httplib

%changelog
* Fri Jun 02 2023 Petr Menšík <pemensik@redhat.com> - 0.12.5-1
- Move back to old-style spec
