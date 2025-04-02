// SPDX-FileCopyrightText: 2025 Pairinteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include <spdlog/details/null_mutex.h>
#include <spdlog/sinks/base_sink.h>
#include <spdlog/spdlog.h>
#include <tbb/concurrent_queue.h>
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
    tbb::concurrent_queue<LogEntry> log_queue;

    class QueueSink : public spdlog::sinks::base_sink<spdlog::details::null_mutex> {
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
