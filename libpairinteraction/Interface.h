#pragma once

#include <string>

int thread_ctrl(int num_threads = -1);

int compute(std::string const &config_name, std::string const &output_name,
            std::string const &wendpoint);

int mainMatrixElement(std::string const &element, std::string const &row,
                      std::string const &col, int power);
