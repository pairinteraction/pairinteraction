%global user            gabime
%global debug_package   %{nil}

Name:           spdlog
Version:        1.7.0
Release:        1%{?dist}
Summary:        Super fast C++ logging library
Group:          Development/Libraries
License:        MIT
URL:            https://github.com/%{user}/%{name}/
Source0:        https://github.com/%{user}/%{name}/archive/v%{version}.tar.gz
BuildRequires:  gcc-c++
BuildRequires:  cmake3
BuildRequires:  fmt-devel

%description
This is a packaged version of the gabime/spdlog header-only C++
logging library available at Github.

%package devel
Summary:        Development files for %{name}
Group:          Development/Libraries
Provides:       %{name}-static = %{version}-%{release}
Provides:       %{name} = %{version}-%{release}
Requires:       libstdc++-devel
Requires:       fmt-devel

%description devel
The %{name}-devel package contains C++ header files for developing
applications that use %{name}.

%prep
%autosetup
mkdir build

%build
pushd build
%{cmake3} .. \
    -DCMAKE_BUILD_TYPE:STRING=RelWithDebInfo \
    -DSPDLOG_FMT_EXTERNAL:BOOL=ON \
    -DSPDLOG_BUILD_EXAMPLE:BOOL=OFF
%make_build

%install
%make_install -C build

%check
pushd build
%{ctest3} -VV %{?_smp_mflags}
popd

%files devel
%{_includedir}/spdlog/
%{_libdir}/libspdlog.*
%{_libdir}/cmake/spdlog/
%{_libdir}/pkgconfig/spdlog.pc

%changelog
* Tue Apr 16 2024 Henri Menke <henri@henrimenke.de> - 1.7.0-1
- Update to 1.7.0

* Sun Sep 04 2016 Daniel Kopecek <dkopecek@redhat.com> - 0.10.0-1
- Update to 0.10.0

* Thu Apr 30 2015 Daniel Kopecek <dkopecek@redhat.com> - 0-4.20150410git211ce99
- don't build the base package
- remove a dot from the release tag
- corrected -devel subpackage description

* Mon Apr 20 2015 Daniel Kopecek <dkopecek@redhat.com> - 0-3.20150410git.211ce99
- use the -p option when copying the header files

* Tue Apr 14 2015 Daniel Kopecek <dkopecek@redhat.com> - 0-2.20150410git.211ce99
- don't build the debuginfo subpackage
- require libstdc++-devel
- don't generate a distribution specific pkg-config file

* Fri Apr 10 2015 Daniel Kopecek <dkopecek@redhat.com> - 0-1.20150410git.211ce99
- Initial package
