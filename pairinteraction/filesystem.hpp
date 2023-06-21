#ifndef FILESYSTEM_H
#define FILESYSTEM_H

// clang-format off
#if __has_include(<filesystem>)
#    include <filesystem>
#elif __has_include(<experimental/filesystem>)
#    include <experimental/filesystem>
#else
#    error "Compiler does not support <filesystem> or <experimental/filesystem>"
#endif
// clang-format on

#include <cereal/types/string.hpp>

#include <memory>

#include <cstdio>
#include <cstdlib>
#include <cstring>

namespace fs {

#if __has_include(<filesystem>)
using namespace std::filesystem;
#elif __has_include(<experimental/filesystem>)
using namespace std::experimental::filesystem;
#endif

/// \brief Generate temporary directory
///
/// \returns path to a temporary directory
/// \warning Not thread-safe!
inline path create_temp_directory() {
#if defined(__unix__) || defined(__linux__) || defined(__APPLE__)
    path const p = temp_directory_path() / "tmpXXXXXX";
    std::unique_ptr<char, decltype(&std::free)> tmpl{strdup(p.c_str()), &std::free};
    char const *const dirname = mkdtemp(tmpl.get());
    if (dirname == nullptr) {
        throw std::runtime_error(std::strerror(errno));
    }
#elif defined(_WIN32)
    char const *const dirname = std::tmpnam(nullptr);
    if (dirname == nullptr) {
        throw std::runtime_error("Could not create temporary directory");
    }
    fs::create_directory(dirname);
#else
#error "Compiler does support not mkdtemp()"
#endif
    return dirname;
}

template <class Archive>
void CEREAL_LOAD_MINIMAL_FUNCTION_NAME(Archive const & /* unused */, path &p,
                                       std::string const &s) {
    p = path(s);
}

template <class Archive>
std::string CEREAL_SAVE_MINIMAL_FUNCTION_NAME(Archive const & /* unused */, path const &p) {
    return p.string();
}

} // namespace fs

#endif // FILESYSTEM_H
