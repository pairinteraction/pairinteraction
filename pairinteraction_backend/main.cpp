#define DOCTEST_CONFIG_IMPLEMENTATION_IN_DLL
#include <doctest/doctest.h>

#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/spdlog.h>

int main(int argc, char **argv) {

    // Configure a logger for the tests
    auto doctest_logger = spdlog::stdout_color_mt("doctest");
    doctest_logger->set_pattern(
        "\033[36m[doctest]\033[0m [%Y-%m-%d %H:%M:%S.%e %t] [%^%l%$] [%s:%#] %v");

    // Run the tests
    doctest::Context ctx;
    ctx.setOption("abort-after", 5);
    ctx.setOption("no-run", 0);
    ctx.applyCommandLine(argc, argv);
    ctx.setOption("no-breaks", true);
    return ctx.run();
}
