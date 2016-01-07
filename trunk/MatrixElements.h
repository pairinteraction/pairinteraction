#ifndef MATRIXELEMENTS_H
#define MATRIXELEMENTS_H

#include "dtypes.h"
#include "Basisnames.h"
#include "DipoleMatrix.hpp"

#include <string>
#include <unordered_map>
#include <gsl/gsl_sf_coupling.h>
#include <memory>

bool selectionRulesDipole(StateOne state1, StateOne state2, int q);
bool selectionRulesMomentum(StateOne state1, StateOne state2, int q);

class MatrixElements {
public:
    MatrixElements(std::string species, int k);
    void precalculate(std::shared_ptr<const BasisnamesOne> basis_one, bool exist_0, bool exist_p, bool exist_m);
    void precalculate(std::shared_ptr<const BasisnamesOne> basis_one, bool exist_d_0, bool exist_d_p, bool exist_d_m, bool exist_m_0, bool exist_m_p, bool exist_m_m);
    real_t getDipole(StateOne state_row, StateOne state_col);
    real_t getMomentum(StateOne state_row, StateOne state_col);
private:
    std::string species;
    int k;
    std::unordered_map<StateTwo,real_t> element_nlj;
    std::unordered_map<StateTwo,real_t> element_nlj_m;
    std::unordered_map<StateTwo,real_t> element_lj_s;
    std::unordered_map<StateTwo,real_t> element_lj_l;
    std::unordered_map<StateTwo,real_t> element_jm;

    real_t muB; // TODO define them in constants.h
    real_t gS;
    real_t gL;
};

#endif
