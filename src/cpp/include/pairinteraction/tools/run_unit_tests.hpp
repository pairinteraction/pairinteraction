// SPDX-FileCopyrightText: 2024 Pairinteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include <filesystem>

namespace pairinteraction {
int run_unit_tests(int argc = 0, char **argv = {}, bool download_missing = false,
                   bool use_cache = true, std::filesystem::path database_dir = "");
} // namespace pairinteraction
