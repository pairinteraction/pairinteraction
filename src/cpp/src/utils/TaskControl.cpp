// SPDX-FileCopyrightText: 2024 PairInteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "pairinteraction/utils/TaskControl.hpp"

#include <atomic>
#include <mutex>

namespace pairinteraction {
namespace {

std::atomic_bool &task_abort_requested() {
    static std::atomic_bool value{false};
    return value;
}

std::mutex &task_info_mutex() {
    static std::mutex value;
    return value;
}

std::string &task_info() {
    static std::string value;
    return value;
}

std::atomic<std::size_t> &progress_count() {
    static std::atomic<std::size_t> value{0};
    return value;
}

} // namespace

TaskAbortedError::TaskAbortedError() : std::runtime_error("Task aborted.") {}

void request_task_abort() noexcept { task_abort_requested().store(true); }

void reset_task_status() noexcept {
    task_abort_requested().store(false);
    progress_count().store(0);
    std::scoped_lock lock(task_info_mutex());
    task_info().clear();
}

std::string get_task_info() {
    std::scoped_lock lock(task_info_mutex());
    return task_info();
}

void set_task_status(std::string_view status_message, bool increase_progress_count) {
    if (!status_message.empty()) {
        std::scoped_lock lock(task_info_mutex());
        task_info().assign(status_message);
    }
    if (increase_progress_count) {
        progress_count().fetch_add(1);
    }
    if (task_abort_requested().load()) {
        throw TaskAbortedError();
    }
}

std::size_t get_progress_count() noexcept { return progress_count().load(); }

} // namespace pairinteraction
