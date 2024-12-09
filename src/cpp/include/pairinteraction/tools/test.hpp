#pragma once

#include <filesystem>

namespace pairinteraction {
int test(int argc = 0, char **argv = {}, bool download_missing = false,
         bool wigner_in_memory = true, std::filesystem::path database_dir = "");
} // namespace pairinteraction
