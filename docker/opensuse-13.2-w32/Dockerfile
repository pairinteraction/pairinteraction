FROM opensuse:13.2
MAINTAINER Pairinteraction https://pairinteraction.github.io
RUN zypper -n --gpg-auto-import-keys addrepo http://download.opensuse.org/repositories/windows:/mingw:/win32/openSUSE_13.2/ windows_mingw_win32 && \
  zypper -n --gpg-auto-import-keys refresh && \
  zypper -n --gpg-auto-import-keys install --no-recommends \
    patterns-openSUSE-devel_C_C++ \
    git cmake curl p7zip \
    eigen3-devel \
    python3-qt5-devel \
    mingw32-cross-gcc-c++ \
    mingw32-boost-devel mingw32-libboost_filesystem mingw32-libboost_program_options mingw32-libboost_system \
    mingw32-libgsl-devel mingw32-libgsl0 \
    sqlite3 mingw32-sqlite-devel mingw32-libsqlite3-0 \
    mingw32-cross-nsis \
    glibc-locale && \
  zypper clean && \
  ln -s /usr/include/eigen3 /usr/i686-w64-mingw32/sys-root/mingw/include/eigen3
