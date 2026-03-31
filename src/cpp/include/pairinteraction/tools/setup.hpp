// SPDX-FileCopyrightText: 2024 PairInteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include <filesystem>

namespace pairinteraction {
void setup(std::filesystem::path ca_bundle_path = {});
} // namespace pairinteraction
