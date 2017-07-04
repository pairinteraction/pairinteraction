#pragma once

#include <string>

int thread_ctrl(int num_threads = -1);

int compute(std::wstring const& config_name, std::wstring const& output_name);

int mainMatrixElement(std::string const& element, std::string const& row, std::string const& col, int power);
