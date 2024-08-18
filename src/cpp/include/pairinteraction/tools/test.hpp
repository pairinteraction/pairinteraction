#pragma once

#include <filesystem>

namespace pairinteraction {
int test(int argc = 0, char **argv = {}, bool download_missing = false,
         std::filesystem::path databasedir = "");
} // namespace pairinteraction
