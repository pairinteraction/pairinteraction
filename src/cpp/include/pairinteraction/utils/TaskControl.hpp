// SPDX-FileCopyrightText: 2024 PairInteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include <stdexcept>
#include <string>
#include <string_view>

namespace pairinteraction {

class TaskAbortedError : public std::runtime_error {
public:
    TaskAbortedError();
};

void request_task_abort() noexcept;

void clear_task_abort() noexcept;

std::string get_task_status();

void task_checkpoint(std::string_view status_message);

} // namespace pairinteraction
