#define DOCTEST_CONFIG_IMPLEMENT

#include "test.hpp"
#include "utils/paths.hpp"

#include <filesystem>
#include <httplib.h>
#include <spdlog/sinks/basic_file_sink.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/spdlog.h>

std::once_flag flag;

int test(int argc, char **argv) {

    // Configure a logger for the tests
    std::filesystem::path logdir = paths::get_pairinteraction_cache_directory() / "logs";

    if (!std::filesystem::exists(logdir)) {
        std::filesystem::create_directories(logdir);
    } else if (!std::filesystem::is_directory(logdir)) {
        throw std::runtime_error("Log path is not a directory.");
    }

    std::filesystem::path logfile = logdir / "test.log";

    std::call_once(flag, [&logfile] {
        auto console_sink = std::make_shared<spdlog::sinks::stdout_color_sink_mt>();
        auto file_sink =
            std::make_shared<spdlog::sinks::basic_file_sink_mt>(logfile.string(), true);
        auto doctest_logger =
            std::make_shared<spdlog::logger>(spdlog::logger("doctest", {console_sink, file_sink}));
        doctest_logger->set_pattern("[doctest] [%Y-%m-%d %H:%M:%S.%e %t] [%^%l%$] [%s:%#] %v");
        spdlog::register_logger(doctest_logger);
    });

    // Setup the tests
    doctest::Context ctx;
    ctx.setOption("abort-after", 5);
    ctx.setOption("no-run", 0);
    ctx.applyCommandLine(argc, argv);
    ctx.setOption("no-breaks", true);

    // Run the tests
    int exitcode = ctx.run();

    if (exitcode != 0) {
        httplib::Client client("https://www.google.com");
        auto res = client.Head("/");
        if (!res) {
            SPDLOG_ERROR("Test failed. Please check your internet connection. An internet "
                         "connection is required to download databases of atomic states and matrix "
                         "elements if they are not available locally.");
        } else {
            SPDLOG_ERROR("Tests failed. Consider creating an issue on "
                         "https://github.com/pairinteraction/pairinteraction/issues, attaching the "
                         "log file {}.",
                         logfile.string());
        }
    }

    return exitcode;
};
