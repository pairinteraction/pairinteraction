// SPDX-FileCopyrightText: 2024 Pairinteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "pairinteraction/tools/setup.hpp"

#include "pairinteraction/utils/paths.hpp"

#include <cstdlib>
#include <filesystem>
#include <mutex>
#include <spdlog/async.h>
#include <spdlog/sinks/rotating_file_sink.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/spdlog.h>
#include <vector>

namespace pairinteraction {
void setup() {

    // Configure a logger
    std::filesystem::path logdir = paths::get_cache_directory() / "logs";

    if (!std::filesystem::exists(logdir)) {
        std::filesystem::create_directories(logdir);
    } else if (!std::filesystem::is_directory(logdir)) {
        throw std::runtime_error("Log path is not a directory.");
    }

    std::filesystem::path logfile = logdir / "cpp.log";

    static std::once_flag flag_default_logger;
    std::call_once(flag_default_logger, [&logfile] {
        spdlog::init_thread_pool(8192, 1);
        auto stdout_sink =
            std::make_shared<spdlog::sinks::stdout_color_sink_mt>(spdlog::color_mode::always);
        auto file_sink = std::make_shared<spdlog::sinks::rotating_file_sink_mt>(logfile.string(),
                                                                                1048576 * 5, 10);
        std::vector<spdlog::sink_ptr> sinks{stdout_sink, file_sink};
        auto logger = std::make_shared<spdlog::async_logger>("logger", sinks.begin(), sinks.end(),
                                                             spdlog::thread_pool(),
                                                             spdlog::async_overflow_policy::block);
        logger->set_pattern("[%Y-%m-%d %H:%M:%S.%e %t] [%^%l%$] [%s:%#] %v");
        spdlog::set_default_logger(logger);

        const char *log_level = std::getenv("SPDLOG_LEVEL");
        if (log_level != nullptr) {
            spdlog::set_level(spdlog::level::from_str(log_level));
        }
    });
}
} // namespace pairinteraction
