#pragma once

#include <filesystem>

namespace pairinteraction {
int run_unit_tests(int argc = 0, char **argv = {}, bool download_missing = false,
                   bool use_cache = true, std::filesystem::path database_dir = "");
} // namespace pairinteraction
