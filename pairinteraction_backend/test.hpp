#pragma once

#include <filesystem>
#include <optional>

int test(int argc = 0, char **argv = {},
         std::optional<bool> download_missing = std::optional<bool>(),
         std::filesystem::path databasedir = "");
