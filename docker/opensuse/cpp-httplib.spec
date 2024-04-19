#
# spec file for package cpp-httplib
#
# Copyright (c) 2023 SUSE LLC
#
# All modifications and additions to the file contributed by third parties
# remain the property of their copyright owners, unless otherwise agreed
# upon. The license for this file, and modifications and additions to the
# file, is the same license as for the pristine package itself (unless the
# license for the pristine package is not an Open Source License, in which
# case the license is the MIT License). An "Open Source License" is a
# license that conforms to the Open Source Definition (Version 1.9)
# published by the Open Source Initiative.

# Please submit bugfixes or comments via https://bugs.opensuse.org/
#


%define sover 0.12
%define libver 0_12
Name:           cpp-httplib
Version:        0.12.5
Release:        0
Summary:        A C++11 HTTP/HTTPS library
License:        MIT
URL:            https://github.com/yhirose/cpp-httplib
Source0:        https://codeload.github.com/yhirose/cpp-httplib/tar.gz/refs/tags/v%{version}#/%{name}-%{version}.tar.gz
BuildRequires:  gcc-c++
BuildRequires:  meson >= 0.47.0
BuildRequires:  pkgconfig(gtest)
BuildRequires:  pkgconfig(libbrotlidec)
BuildRequires:  pkgconfig(libbrotlienc)
BuildRequires:  pkgconfig(openssl) >= 1.1.1

%package -n lib%{name}%{libver}
Summary:        A C++11 HTTP/HTTPS library

%package devel
Summary:        A C++11 HTTP/HTTPS library
Requires:       lib%{name}%{libver} = %{version}
Conflicts:      cpp-httplib-headers-devel

%description
This is a multi-threaded HTTP library with blocking I/O. There is no
support for non-blocking mode.

%description -n lib%{name}%{libver}
This is a multi-threaded HTTP library with blocking I/O. There is no
support for non-blocking mode.

%description devel
This is a multi-threaded HTTP library with blocking I/O. There is no
support for non-blocking mode.

It features built-in mappings, static file server, pre-routing and
post-routing handlers, and support for binding sockets to multiple
interfaces and any available port.

%prep
%setup -q
chmod -x example/uploader.sh

%build
%meson -Dcpp-httplib_compile=true -Dcpp-httplib_test=true \
       --buildtype=release
%meson_build

%install
%meson_install

%check
# OBS and chroot build environments does not provide internet
# connectivity, skip online tests to avoid failures
export GTEST_FILTER='-*.*_Online'
%meson_test

%post   -n lib%{name}%{libver} -p /sbin/ldconfig
%postun -n lib%{name}%{libver} -p /sbin/ldconfig

%files -n lib%{name}%{libver}
%{_libdir}/lib%{name}.so.%{sover}
%{_libdir}/lib%{name}.so.%{version}
%license LICENSE

%files devel
%{_libdir}/lib%{name}.so
%{_includedir}/httplib.h
%{_libdir}/pkgconfig/%{name}.pc
%doc README.md example

%changelog
