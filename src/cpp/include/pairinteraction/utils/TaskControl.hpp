// SPDX-FileCopyrightText: 2024 PairInteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include <cstddef>
#include <stdexcept>
#include <string>
#include <string_view>

namespace pairinteraction {

class TaskAbortedError : public std::runtime_error {
public:
    TaskAbortedError();
};

void request_task_abort() noexcept;

void reset_task_status() noexcept;

std::string get_task_info();

void set_task_status(std::string_view status_message, bool increase_progress_count = false);

std::size_t get_progress_count() noexcept;

} // namespace pairinteraction
