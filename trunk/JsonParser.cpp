#include "JsonParser.hpp"
#include <string>
#include <fstream>
#include <iostream>
#include "json/json.h"


JsonParser::JsonParser(std::string filename)
  : filename(filename) {}

int JsonParser::parse()
{
  std::ifstream conffile(filename);
  if ( !reader.parse(conffile, datadict) )
  {
    std::cerr << reader.getFormattedErrorMessages();
    return 1;
  }
  
  element            = datadict["element"].asString();
  n1                 = datadict["n1"].asInt();
  n2                 = datadict["n2"].asInt();
  l1                 = datadict["l1"].asInt();
  l2                 = datadict["l2"].asInt();
  j1                 = datadict["j1"].asDouble();
  j2                 = datadict["j2"].asDouble();
  m1                 = datadict["m1"].asDouble();
  m2                 = datadict["m2"].asDouble();
  symmetry           = datadict["symmetry"].asString();
  angle              = datadict["angle"].asDouble();
  enable_dipdip      = datadict["enable_dipdip"].asBool();
  enable_dipquad     = datadict["enable_dipquad"].asBool();
  enable_quadquad    = datadict["enable_quadquad"].asBool();
  efield_strength    = datadict["efield_strength"].asDouble();
  efield_increasing  = datadict["efield_increasing"].asBool();
  bfield_strength    = datadict["bfield_strength"].asDouble();
  bfield_increasing  = datadict["bfield_increasing"].asBool();
  delta_n            = datadict["delta_n"].asInt();
  delta_l            = datadict["delta_l"].asInt();
  delta_m            = datadict["delta_m"].asInt();
  delta_energy       = datadict["delta_energy"].asDouble();
  preserve_M         = datadict["preserve_M"].asBool();
  preserve_submatrix = datadict["preserve_submatrix"].asBool();
  preserve_parityL   = datadict["preserve_parityL"].asBool();
  distance_min       = datadict["distance_min"].asDouble();
  distance_max       = datadict["distance_max"].asDouble();
  distance_steps     = datadict["distance_steps"].asInt();
  energyrange        = datadict["energyrange"].asDouble();

  return 0;
}
