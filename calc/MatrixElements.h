#ifndef MATRIXELEMENTS_H
#define MATRIXELEMENTS_H

#include "dtypes.h"
#include "Basisnames.h"
#include "Numerov.hpp"

#include <string>
#include <unordered_map>
//#include <gsl/gsl_sf_coupling.h>
#include <memory>
#include <sstream>
#include <wignerSymbols.h>

bool selectionRulesDipole(StateOne state1, StateOne state2, int q);
bool selectionRulesQuadrupole(StateOne state1, StateOne state2, int q);
bool selectionRulesMomentum(StateOne state1, StateOne state2, int q);

size_t findidx(std::vector<real_t> x, real_t d);

class MatrixElements { // TODO ein buffer am Programmstart
public:
    MatrixElements(std::string species, int k, std::string dbname);
    void precalculate_momentum(std::shared_ptr<const BasisnamesOne> basis_one, bool exist_0, bool exist_p, bool exist_m);
    void precalculate_dipole(std::shared_ptr<const BasisnamesOne> basis_one, bool exist_0, bool exist_p, bool exist_m);
    void precalculate_quadrupole(std::shared_ptr<const BasisnamesOne> basis_one, bool exist_0, bool exist_p, bool exist_m, bool exist_pp, bool exist_mm);
    real_t getDipole(StateOne state_row, StateOne state_col);
    real_t getQuadrupole(StateOne state_row, StateOne state_col);
    real_t getMomentum(StateOne state_row, StateOne state_col);
private:
    void precalculate(std::shared_ptr<const BasisnamesOne> basis_one, bool exist_d_0, bool exist_d_p, bool exist_d_m, bool exist_q_0, bool exist_q_p, bool exist_q_m, bool exist_q_pp, bool exist_q_mm, bool exist_m_0, bool exist_m_p, bool exist_m_m);
    real_t calcRadialElement(std::string species, int n1, int l1, real_t j1, int power, int n2, int l2, real_t j2);
    std::string species;
    int k;
    std::string dbname;
    std::unordered_map<StateTwo,real_t> element_nlj_k;
    std::unordered_map<StateTwo,real_t> element_nlj_0;
    std::unordered_map<StateTwo,real_t> element_nlj_m;
    std::unordered_map<StateTwo,real_t> element_lj_s;
    std::unordered_map<StateTwo,real_t> element_lj_l;
    std::unordered_map<StateTwo,real_t> element_jm;

    real_t muB; // TODO define them in constants.h
    real_t gS;
    real_t gL;
};

#endif
