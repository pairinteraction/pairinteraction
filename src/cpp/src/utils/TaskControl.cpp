// SPDX-FileCopyrightText: 2024 PairInteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "pairinteraction/utils/TaskControl.hpp"

#include <atomic>
#include <mutex>

namespace pairinteraction {
std::atomic_bool &task_abort_requested() {
    static std::atomic_bool value = false;
    return value;
}

std::mutex &task_status_mutex() {
    static std::mutex value;
    return value;
}

std::string &task_status() {
    static std::string value;
    return value;
}

TaskAbortedError::TaskAbortedError() : std::runtime_error("Task aborted.") {}

void request_task_abort() noexcept {
    task_abort_requested().store(true, std::memory_order_relaxed);
}

void clear_task_abort() noexcept {
    task_abort_requested().store(false, std::memory_order_relaxed);
    std::scoped_lock lock(task_status_mutex());
    task_status().clear();
}

std::string get_task_status() {
    std::scoped_lock lock(task_status_mutex());
    return task_status();
}

void task_checkpoint(std::string_view status_message) {
    if (!status_message.empty()) {
        std::scoped_lock lock(task_status_mutex());
        task_status().assign(status_message);
    }
    if (task_abort_requested().load(std::memory_order_relaxed)) {
        throw TaskAbortedError();
    }
}

} // namespace pairinteraction
