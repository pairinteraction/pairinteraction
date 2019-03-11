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
#include "Interface.h"

#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <iostream>
#include <memory>

int main(int argc, char **argv) {
    namespace po = boost::program_options;

    po::options_description desc("Usage");
    desc.add_options()("help,?", "produce this help message")(
        "config,c", po::value<std::string>()->required(), "Path to config JSON file")(
        "output,o", po::value<std::string>()->required(), "Path to cache JSON file");

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);

    if (vm.count("help") != 0u) {
        std::cout << desc << std::endl;
        return 0;
    }

    try {
        po::notify(vm);
    } catch (po::required_option &e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    return compute(vm["config"].as<std::string>(), vm["output"].as<std::string>());
}
