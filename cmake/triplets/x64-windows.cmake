# SPDX-FileCopyrightText: 2024 PairInteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

set(VCPKG_TARGET_ARCHITECTURE x64)
set(VCPKG_CRT_LINKAGE dynamic)
set(VCPKG_LIBRARY_LINKAGE dynamic)
set(VCPKG_BUILD_TYPE release)

# Link OpenSSL statically. Python ships its own OpenSSL DLLs under the same names (libssl-3-x64.dll,
# libcrypto-3-x64.dll). If Python has already loaded them, e.g. because the ssl module was imported before
# pairinteraction, Windows reuses them for our backend and loading fails if their version is older.
if(PORT STREQUAL "openssl")
  set(VCPKG_LIBRARY_LINKAGE static)
endif()
