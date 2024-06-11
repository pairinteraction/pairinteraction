#pragma once

#include <filesystem>

int test(int argc = 0, char **argv = {}, bool download_missing = false,
         std::filesystem::path databasedir = "");
