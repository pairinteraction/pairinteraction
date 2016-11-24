FROM opensuse:13.2
MAINTAINER Henri Menke <pairinteraction@henrimenke.de>

ENV TERM ansi

RUN zypper -n --gpg-auto-import-keys addrepo http://download.opensuse.org/repositories/windows:/mingw:/win32/openSUSE_13.2/ windows_mingw_win32 && \
  zypper -n --gpg-auto-import-keys refresh && \
  zypper -n --gpg-auto-import-keys install --no-recommends \
    patterns-openSUSE-devel_C_C++ \
    git cmake curl p7zip \
    eigen3-devel \
    python3-qt5-devel \
    mingw32-cross-gcc-c++ \
    mingw32-boost-devel mingw32-libboost_filesystem mingw32-libboost_system \
    mingw32-libgsl-devel \
    sqlite3 mingw32-sqlite-devel \
    mingw32-cross-nsis && \
  zypper clean && \
  ln -s /usr/include/eigen3 /usr/i686-w64-mingw32/sys-root/mingw/include/eigen3

RUN zypper --non-interactive --gpg-auto-import-keys install curl
RUN zypper --non-interactive --gpg-auto-import-keys install sudo
RUN zypper --non-interactive --gpg-auto-import-keys install mingw32-cross-nsis-plugin-uac mingw32-cross-nsis-plugin-nsprocess
RUN zypper --non-interactive --gpg-auto-import-keys install mingw32-zlib-devel
RUN zypper --non-interactive --gpg-auto-import-keys install glibc-locale
