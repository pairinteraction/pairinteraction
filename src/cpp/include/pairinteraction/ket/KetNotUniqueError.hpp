// SPDX-FileCopyrightText: 2026 PairInteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include <memory>
#include <stdexcept>
#include <utility>
#include <vector>

namespace pairinteraction {
class KetAtom;

/**
 * @class KetNotUniqueError
 *
 * @brief Exception thrown when a ket is not uniquely specified.
 *
 * The exception carries the candidate kets that match the description so that the caller (e.g. the
 * Python bindings) can format a helpful error message using its own labeling logic.
 */
class KetNotUniqueError : public std::runtime_error {
public:
    explicit KetNotUniqueError(std::vector<std::shared_ptr<const KetAtom>> kets)
        : std::runtime_error("The ket is not uniquely specified."), kets(std::move(kets)) {}

    const std::vector<std::shared_ptr<const KetAtom>> &get_kets() const { return kets; }

private:
    std::vector<std::shared_ptr<const KetAtom>> kets;
};

} // namespace pairinteraction
