#ifndef MATRIXELEMENTS_H
#define MATRIXELEMENTS_H

#include "dtypes.h"
#include "Basisnames.h"
#include "Wavefunction.hpp"

#include <string>
#include <unordered_map>
//#include <gsl/gsl_sf_coupling.h>
#include <memory>
#include <sstream>
#include <wignerSymbols.h>

bool selectionRulesMultipole(StateOne state1, StateOne state2, int kappa);
bool selectionRulesMomentum(StateOne state1, StateOne state2);


bool selectionRulesDipole(StateOne state1, StateOne state2, int q);
bool selectionRulesQuadrupole(StateOne state1, StateOne state2, int q);
bool selectionRulesMomentumOld(StateOne state1, StateOne state2, int q);

class MatrixElements { // TODO ein buffer am Programmstart
public:
    MatrixElements(std::string species, std::string dbname);
    MatrixElements(const Configuration& config, std::string species, std::string dbname);
    void precalculate_momentumOld(std::shared_ptr<const BasisnamesOne> basis_one, bool exist_0, bool exist_p, bool exist_m);
    void precalculateElectricMomentum(std::shared_ptr<const BasisnamesOne> basis_one, int q);
    void precalculate_quadrupole(std::shared_ptr<const BasisnamesOne> basis_one, bool exist_0, bool exist_p, bool exist_m, bool exist_pp, bool exist_mm);
    void precalculate_multipole(std::shared_ptr<const BasisnamesOne> basis_one, int kappa);
    void precalculateMagneticMomentum(std::shared_ptr<const BasisnamesOne> basis_one, int q);
    real_t getMultipole(StateOne state_row, StateOne state_col, int k);
    real_t getMagneticMomentum(StateOne state_row, StateOne state_col);
    real_t getElectricMomentum(StateOne state_row, StateOne state_col);

    real_t getQuadrupole(StateOne state_row, StateOne state_col);
    real_t getMomentumOld(StateOne state_row, StateOne state_col);
private:
    void precalculateOld(std::shared_ptr<const BasisnamesOne> basis_one, bool exist_d_0, bool exist_d_p, bool exist_d_m, bool exist_q_0, bool exist_q_p, bool exist_q_m, bool exist_q_pp, bool exist_q_mm, bool exist_m_0, bool exist_m_p, bool exist_m_m);
    void precalculate(std::shared_ptr<const BasisnamesOne> basis_one, int kappa, int q, bool calcMultipole, bool calcMomentum);
    real_t calcRadialElement(std::string species, int n1, int l1, real_t j1, int power, int n2, int l2, real_t j2);
    std::string method;
    std::string species;
    std::string dbname;
    std::unordered_map<StateTwo,real_t> element_nlj_k;
    std::unordered_map<StateTwo,real_t> element_nlj_0;
    std::unordered_map<StateTwo,real_t> element_nlj_m;
    std::unordered_map<StateTwo,real_t> element_lj_s;
    std::unordered_map<StateTwo,real_t> element_lj_l;
    std::unordered_map<StateTwo,real_t> element_jm;

    std::unordered_map<int,std::unordered_map<StateTwo,real_t>> cache_radial;
    std::unordered_map<int,std::unordered_map<StateTwo,real_t>> cache_angular;
    std::unordered_map<int,std::unordered_map<StateTwo,real_t>> cache_reduced_commutes_s;
    std::unordered_map<int,std::unordered_map<StateTwo,real_t>> cache_reduced_commutes_l;
    std::unordered_map<int,std::unordered_map<StateTwo,real_t>> cache_reduced_multipole;

    real_t muB; // TODO define them in constants.h
    real_t gS;
    real_t gL;
};

#endif
