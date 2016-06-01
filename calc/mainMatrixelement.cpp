#include "MatrixElements.h"
#include "QuantumDefect.hpp"
#include "State.h"
#include "Basisnames.h"

#include <getopt.h>
#include <sstream>
#include <string>
#include <limits>
#include <boost/filesystem.hpp>

// ############################################################################
// ### MAIN LOOP ##############################################################
// ############################################################################

int main(int argc, char **argv) {
    std::cout << std::unitbuf;

    std::string str_row = "";
    std::string str_col = "";

    int c;
    opterr = 0;
    while((c = getopt (argc, argv, "r:c:")) != -1) {
        switch (c) {
        case 'r':
            str_row = optarg;
            break;
        case 'c':
            str_col = optarg;
            break;
        default:
            return 1;
        }
    }

    size_t precision = std::numeric_limits<real_t>::digits10 + 1;

    if (str_row != "" and str_col != "") { // Dipole matrix element
        std::stringstream ss_row(str_row);
        std::stringstream ss_col(str_col);

        StateOne state_row, state_col;
        state_row.m = 0;
        state_col.m = 0;
        ss_row >> state_row.n >> state_row.l >> state_row.j >> state_row.m;
        ss_col >> state_col.n >> state_col.l >> state_col.j >> state_col.m;

        real_t dipolematrixelement = 0;
        if (selectionRulesDipole(state_row, state_col, state_row.m-state_col.m)) {
            MatrixElements matrixelement("Rb", 1, "");
            std::vector<StateOne> states({{state_row, state_col}});
            auto basis = std::make_shared<BasisnamesOne>(BasisnamesOne::fromStates(states));
            matrixelement.precalculate_dipole(basis, true, true, true);
            dipolematrixelement = matrixelement.getDipole(state_row, state_col);
        }

        std::cout << std::setprecision(precision) << dipolematrixelement << std::endl;

    } else if (str_row != "" or str_col != "") {
        std::stringstream ss;
        if (str_row != "") { // Energy of row state
            ss = std::stringstream(str_row);
        } else { // Energy of columne state
            ss = std::stringstream(str_col);
        }

        StateOne state;
        ss >> state.n >> state.l >> state.j >> state.m;

        real_t energy = energy_level("Rb", state.n, state.l, state.j);

        std::cout << std::setprecision(precision) << energy << std::endl;
    }

    return 0;
}
