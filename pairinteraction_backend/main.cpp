#include "setup.hpp"
#include "test.hpp"

#ifndef PAIRINTERACTION_LOCAL_DATABASE_DIR
#error PAIRINTERACTION_LOCAL_DATABASE_DIR is not defined
#endif

int main(int argc, char **argv) {
    setup();
    return test(argc, argv, false, PAIRINTERACTION_LOCAL_DATABASE_DIR);
}
