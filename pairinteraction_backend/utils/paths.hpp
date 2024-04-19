#include <cstdlib>
#include <filesystem>

namespace paths {

std::filesystem::path get_home_directory() {
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
    return std::filesystem::path("~");
}

std::filesystem::path get_pairinteraction_cache_directory() {
    std::filesystem::path path;

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

} // namespace paths
