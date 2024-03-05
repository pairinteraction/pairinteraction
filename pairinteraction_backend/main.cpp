#define DOCTEST_CONFIG_IMPLEMENTATION_IN_DLL
#include <doctest/doctest.h>

int main(int argc, char **argv) {

    // Run the tests
    doctest::Context ctx;
    ctx.setOption("abort-after", 5);
    ctx.setOption("no-run", 0);
    ctx.applyCommandLine(argc, argv);
    ctx.setOption("no-breaks", true);
    return ctx.run();
}
