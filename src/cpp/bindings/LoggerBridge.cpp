// SPDX-FileCopyrightText: 2025 Pairinteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "./LoggerBridge.hpp"

#include "pairinteraction/utils/paths.hpp"

#include <filesystem>
#include <fmt/core.h>
#include <memory>
#include <spdlog/async.h>
#include <spdlog/sinks/base_sink.h>
#include <spdlog/sinks/rotating_file_sink.h>
#include <spdlog/spdlog.h>
#include <string>
#include <vector>

using namespace pairinteraction;

LoggerBridge::QueueSink::QueueSink(LoggerBridge *parent) : parent(parent) {}

void LoggerBridge::QueueSink::sink_it_(const spdlog::details::log_msg &msg) {
    spdlog::memory_buf_t buf;
    this->formatter_->format(msg, buf);
    std::string text = fmt::to_string(buf);
    LogEntry entry;
    switch (msg.level) {
    case spdlog::level::trace:
        entry.level = 0;
        break;
    case spdlog::level::debug:
        entry.level = 10;
        break;
    case spdlog::level::info:
        entry.level = 20;
        break;
    case spdlog::level::warn:
        entry.level = 30;
        break;
    case spdlog::level::err:
        entry.level = 40;
        break;
    case spdlog::level::critical:
        entry.level = 50;
        break;
    default:
        break;
    }
    entry.message = text;
    parent->log_queue.push(std::move(entry));
}

void LoggerBridge::QueueSink::flush_() {}

LoggerBridge::LoggerBridge() {
    std::filesystem::path logdir = paths::get_cache_directory() / "logs";
    if (!std::filesystem::exists(logdir)) {
        std::filesystem::create_directories(logdir);
    } else if (!std::filesystem::is_directory(logdir)) {
        throw std::runtime_error("Log path is not a directory.");
    }
    std::filesystem::path logfile = logdir / "cpp.log";

    spdlog::init_thread_pool(8192, 1);

    auto queue_sink = std::make_shared<QueueSink>(this);
    queue_sink->set_pattern("[%Y-%m-%d %H:%M:%S.%e %t] [%s:%#] %v");

    auto file_sink =
        std::make_shared<spdlog::sinks::rotating_file_sink_mt>(logfile.string(), 1048576 * 5, 10);
    file_sink->set_pattern("[%Y-%m-%d %H:%M:%S.%e %t] [%^%l%$] [%s:%#] %v");

    std::vector<spdlog::sink_ptr> sinks{queue_sink, file_sink};
    auto logger = std::make_shared<spdlog::async_logger>("logger", sinks.begin(), sinks.end(),
                                                         spdlog::thread_pool(),
                                                         spdlog::async_overflow_policy::block);
    logger->set_level(spdlog::level::trace);
    spdlog::set_default_logger(logger);
}

LoggerBridge::~LoggerBridge() {
    spdlog::shutdown();
    logger.reset();
}

std::vector<LoggerBridge::LogEntry> LoggerBridge::get_pending_logs() {
    std::vector<LogEntry> entries;
    LogEntry entry;
    while (log_queue.try_pop(entry)) {
        entries.push_back(std::move(entry));
    }
    return entries;
}
