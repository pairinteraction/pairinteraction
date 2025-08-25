// SPDX-FileCopyrightText: 2025 Pairinteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include <deque>
#include <mutex>
#include <spdlog/sinks/base_sink.h>
#include <spdlog/spdlog.h>
#include <string>
#include <vector>

class LoggerBridge {
public:
    struct LogEntry {
        int level = 20;
        std::string message;
    };

    LoggerBridge();
    ~LoggerBridge();

    std::vector<LogEntry> get_pending_logs();

private:
    std::deque<LogEntry> log_queue;
    std::mutex log_queue_mtx;

    class QueueSink : public spdlog::sinks::base_sink<std::mutex> {
    public:
        explicit QueueSink(LoggerBridge *parent);

    protected:
        void sink_it_(const spdlog::details::log_msg &msg) override;
        void flush_() override;

    private:
        LoggerBridge *parent;
    };

    std::shared_ptr<spdlog::logger> logger;
};
