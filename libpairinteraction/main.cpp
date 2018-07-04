/*
 * Copyright (c) 2016 Sebastian Weber, Henri Menke. All rights reserved.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
#include "Interface.h"

#include <iostream>
#include <memory>
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>

int main(int argc, char **argv) {
    namespace po = boost::program_options;

    po::options_description desc("Usage");
    desc.add_options()
            ("help,?", "produce this help message")
            ("config,c", po::value<std::string>()->required(),"Path to config JSON file")
            ("output,o", po::value<std::string>()->required(),"Path to cache JSON file")
            ;

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);

    if ( vm.count("help") != 0u )
    {
        std::cout << desc << std::endl;
        return 0;
    }

    try
    {
        po::notify(vm);
    }
    catch (po::required_option& e)
    {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    return compute(vm["config"].as<std::string>(), vm["output"].as<std::string>(), "use_cout");
}
