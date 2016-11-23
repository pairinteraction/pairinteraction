FROM opensuse:13.2
MAINTAINER Henri Menke <pairinteraction@henrimenke.de>
RUN zypper -n --gpg-auto-import-keys addrepo http://download.opensuse.org/repositories/windows:/mingw:/win64/openSUSE_13.2/ windows_mingw_win64 && \
  zypper -n --gpg-auto-import-keys addrepo http://download.opensuse.org/repositories/windows:/mingw:/win32/openSUSE_13.2/ windows_mingw_win32 && \
  zypper -n --gpg-auto-import-keys refresh && \
  zypper -n --gpg-auto-import-keys install --no-recommends \
    patterns-openSUSE-devel_C_C++ \
    git cmake curl p7zip \
    eigen3-devel \
    python3-qt5-devel \
    mingw64-cross-gcc-c++ \
    mingw64-boost-devel mingw64-libboost_filesystem mingw64-libboost_system \
    mingw64-libgsl-devel \
    sqlite3 mingw64-sqlite-devel
    mingw32-cross-nsis && \
  ln -s /usr/include/eigen3 /usr/x86_64-w64-mingw32/sys-root/mingw/include/eigen3
