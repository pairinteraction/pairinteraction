/*
 * Copyright (c) 2016 Sebastian Weber, Henri Menke. All rights reserved.
 *
 * This file is part of the pairinteraction library.
 *
 * The pairinteraction library is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * The pairinteraction library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with the pairinteraction library. If not, see <http://www.gnu.org/licenses/>.
 */
#include "Interface.hpp"

#include <iostream>
#include <string>

static void print_usage(std::ostream &os, int status) {
    os << "Usage:\n"
          "  -? [ --help ]         produce this help message\n"
          "  -c [ --config ] arg   Path to config JSON file\n"
          "  -o [ --output ] arg   Path to cache JSON file\n"
          "  -t [ --type ] arg     Type of calculation: arg = real|complex\n";
    std::exit(status);
}

int main(int argc, char *argv[]) {
    if (argc == 1) {
        print_usage(std::cout, EXIT_SUCCESS);
    }

    std::string config, output, type;

    int optind = 1;
    while (optind < argc) {
        std::string opt = argv[optind];
        if (opt == "-?" || opt == "--help") {
            print_usage(std::cout, EXIT_SUCCESS);
        } else if (opt == "-c" || opt == "--config") {
            ++optind;
            if (!(optind < argc)) {
                std::cerr << "Option " << opt << " requires an argument\n";
                std::exit(EXIT_FAILURE);
            }
            config = argv[optind];
        } else if (opt == "-o" || opt == "--output") {
            ++optind;
            if (!(optind < argc)) {
                std::cerr << "Option " << opt << " requires an argument\n";
                std::exit(EXIT_FAILURE);
            }
            output = argv[optind];
        } else if (opt == "-t" || opt == "--type") {
            ++optind;
            if (!(optind < argc)) {
                std::cerr << "Option " << opt << " requires an argument\n";
                std::exit(EXIT_FAILURE);
            }
            type = argv[optind];
        } else {
            std::cerr << "Unknown option: " << opt << "\n";
            print_usage(std::cerr, EXIT_FAILURE);
        }
        ++optind;
    }

    if (config.empty()) {
        std::cerr << "Option --config is required\n";
        std::exit(EXIT_FAILURE);
    }

    if (output.empty()) {
        std::cerr << "Option --output is required\n";
        std::exit(EXIT_FAILURE);
    }

    if (type.empty()) {
        std::cerr << "Option --type is required\n";
        std::exit(EXIT_FAILURE);
    }

    if (type == "complex") {
        return compute<std::complex<double>>(config, output);
    } else if (type == "real") {
        return compute<double>(config, output);
    } else {
        std::cerr << "Option --type can only be 'real' or 'complex' but is '" << type << "'\n";
        std::exit(EXIT_FAILURE);
    }
}
