// SPDX-FileCopyrightText: 2024 Pairinteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include <cstdlib>
#include <filesystem>

namespace pairinteraction::paths {

inline std::filesystem::path get_home_directory() {
#ifdef _WIN32
    char const *hdrive = getenv("HOMEDRIVE");
    char const *hpath = getenv("HOMEPATH");
    if (hdrive != nullptr && hpath != nullptr) {
        std::filesystem::path path = std::string(hdrive) + hpath;
        return path;
    }
#else
    const char *home = std::getenv("HOME");
    if (home != nullptr) {
        std::filesystem::path path = home;
        return path;
    }
#endif
    return {"~"};
}

inline std::filesystem::path get_cache_directory() {
    std::filesystem::path path;

    char const *pairinteraction_cache_dir = getenv("PAIRINTERACTION_CACHE_DIR");
    if (pairinteraction_cache_dir != nullptr) {
        path = pairinteraction_cache_dir;
        return path;
    }

#ifdef _WIN32
    char const *localappdata = std::getenv("LOCALAPPDATA");
    if (localappdata != nullptr) {
        path = localappdata;
    } else {
        path = get_home_directory() / "AppData" / "Local";
    }
    path /= "pairinteraction";
#elif __APPLE__
    path = get_home_directory() / "Library" / "Caches" / "pairinteraction";
#else
    char const *xdg_cache_home = std::getenv("XDG_CACHE_HOME");
    if (xdg_cache_home != nullptr) {
        path = xdg_cache_home;
    } else {
        path = get_home_directory() / ".cache";
    }
    path /= "pairinteraction";
#endif

    return path;
}

inline std::filesystem::path get_config_directory() {
    std::filesystem::path path;

    char const *pairinteraction_config_dir = getenv("PAIRINTERACTION_CONFIG_DIR");
    if (pairinteraction_config_dir != nullptr) {
        path = pairinteraction_config_dir;
        return path;
    }

#ifdef _WIN32
    char const *appdata = std::getenv("APPDATA");
    if (appdata != nullptr) {
        path = appdata;
    } else {
        path = get_home_directory() / "AppData" / "Roaming";
    }
    path /= "pairinteraction";
#elif __APPLE__
    path = get_home_directory() / "Library" / "Preferences" / "pairinteraction";
#else
    char const *xdg_config_home = std::getenv("XDG_CONFIG_HOME");
    if (xdg_config_home != nullptr) {
        path = xdg_config_home;
    } else {
        path = get_home_directory() / ".config";
    }
    path /= "pairinteraction";
#endif

    return path;
}

} // namespace pairinteraction::paths
