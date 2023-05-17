#pragma once

#include <complex>
#include <string>

template <typename Scalar>
int compute(std::string const &config_name, std::string const &output_name);

#ifndef SWIG
extern template int compute<std::complex<double>>(std::string const &config_name, std::string const &output_name);
extern template int compute<double>(std::string const &config_name, std::string const &output_name);
#endif
