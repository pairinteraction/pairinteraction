#include "MatrixElements.h"
#include "SQLite.hpp"
#include <sstream>
#include <iostream>
#include <string>

bool selectionRulesDipole(StateOne state1, StateOne state2, int q) {
    return (abs(state1.l-state2.l) == 1) && (state1.m == state2.m+q) && (fabs(state1.j-state2.j) <= 1);
}

bool selectionRulesMomentum(StateOne state1, StateOne state2, int q) {
    return (state1.l == state2.l) && (state1.m == state2.m+q) && (fabs(state1.j-state2.j) <= 1) && (state1.n == state2.n);
}

MatrixElements::MatrixElements() {
    k = 1;
    alkali = "Rb";

    muB = 0.5;
    gS = 2.0023192;
    gL = 1;
}

void MatrixElements::precalculate(BasisnamesOne basis_one, bool exist_d_0, bool exist_d_p, bool exist_d_m, bool exist_m_0, bool exist_m_p, bool exist_m_m) {
    //TODO: check if k = 1 in case of exist_m

    // determine elements

    SQLite3 db("matrix_elements.db");

    db.exec("CREATE TEMPORARY TABLE tmp_nlj ("
            "n1 integer, l1 integer, j1 double,"
            "n2 integer, l2 integer, j2 double);");

    db.exec("CREATE TEMPORARY TABLE tmp_lj_s ("
            "l1 integer, j1 double,"
            "l2 integer, j2 double);");


    db.exec("CREATE TEMPORARY TABLE tmp_lj_l ("
            "l1 integer, j1 double,"
            "l2 integer, j2 double);");


    db.exec("CREATE TEMPORARY TABLE tmp_jm ("
            "j1 integer, m1 double,"
            "j2 integer, m2 double);");

    std::stringstream ss;

    db.exec("begin transaction;");

    for (const auto &state_col : basis_one) {
        for (const auto &state_row : basis_one) {
            if (state_row.idx < state_col.idx) {
                continue;
            }

            if ((exist_d_0 && selectionRulesDipole(state_row, state_col, 0)) ||
                    (exist_d_p && selectionRulesDipole(state_row, state_col, 1)) ||
                    (exist_d_m && selectionRulesDipole(state_row, state_col, -1)) ||
                    (exist_m_0 && selectionRulesMomentum(state_row, state_col, 0)) ||
                    (exist_m_p && selectionRulesMomentum(state_row, state_col, 1)) ||
                    (exist_m_m && selectionRulesMomentum(state_row, state_col, -1))) {

                if ((exist_d_0 || exist_d_p || exist_d_m) && (abs(state_row.l-state_col.l) == 1)) {
                    auto result_nlj = element_nlj.insert(std::make_pair<StateTwo,real_t>(
                                                               StateTwo({{state_row.n, state_col.n}}, {{state_row.l, state_col.l}}, {{0,0}}, {{state_row.j, state_col.j}}, {{0,0}}).order(), std::numeric_limits<real_t>::max()));
                    if (result_nlj.second) {
                        ss.str(std::string());
                        ss << "insert into tmp_nlj (n1,l1,j1,n2,l2,j2) values ("
                           << state_row.n << "," << state_row.l << "," << state_row.j << ","
                           << state_col.n << "," << state_col.l << "," << state_col.j << ");";
                        db.exec(ss.str().c_str());
                    }
                }

                if ((exist_m_0 || exist_m_p || exist_m_m) && (state_row.l == state_col.l) && (state_row.n == state_col.n) ) {
                    auto result_lj = element_lj_s.insert(std::make_pair<StateTwo,real_t>(
                                                           StateTwo({{0,0}}, {{state_row.l, state_col.l}}, {{0,0}}, {{state_row.j, state_col.j}}, {{0,0}}).order(), std::numeric_limits<real_t>::max()));
                    if (result_lj.second) {
                        ss.str(std::string());
                        ss << "insert into tmp_lj_l (l1,j1,l2,j2) values ("
                           << state_row.l << "," << state_row.j << ","
                           << state_col.l << "," << state_col.j << ");";
                        db.exec(ss.str().c_str());
                    }
                }

                auto result_lj = element_lj_l.insert(std::make_pair<StateTwo,real_t>(
                                                       StateTwo({{0,0}}, {{state_row.l, state_col.l}}, {{0,0}}, {{state_row.j, state_col.j}}, {{0,0}}).order(), std::numeric_limits<real_t>::max()));
                if (result_lj.second) {
                    ss.str(std::string());
                    ss << "insert into tmp_lj_l (l1,j1,l2,j2) values ("
                       << state_row.l << "," << state_row.j << ","
                       << state_col.l << "," << state_col.j << ");";
                    db.exec(ss.str().c_str());
                }

                auto result_jm = element_jm.insert(std::make_pair<StateTwo,real_t>(
                                                       StateTwo({{0,0}}, {{0, 0}}, {{0,0}}, {{state_row.j, state_col.j}}, {{state_row.m, state_col.m}}).order(), std::numeric_limits<real_t>::max()));
                if (result_jm.second) {
                    ss.str(std::string());
                    ss << "insert into tmp_jm (j1,m1,j2,m2) values ("
                       << state_row.j << "," << state_row.m << ","
                       << state_col.j << "," << state_col.m << ");";
                    db.exec(ss.str().c_str());
                }
            }
        }
    }

    db.exec("end transaction;");

    std::cout << 8.31 << std::endl;

    // load from database
    // TODO

    // calculate missing elements
    if (exist_d_0 || exist_d_p || exist_d_m) {
        for (auto &element : element_nlj) {
            if (element.second == std::numeric_limits<real_t>::max()) {
                int lmax = fmax(element.first.l[0],element.first.l[1]);

                element.second = pow(-1, element.first.l[0]+lmax) * sqrt(lmax) *
                        radial_element(alkali, element.first.n[0], element.first.l[0], element.first.j[0], k, element.first.n[1], element.first.l[1], element.first.j[1]); // TODO <Rb|Cs> not possible, so Rb is enough
                //TODO: query
            }
        }
    }

    if (exist_m_0 || exist_m_p || exist_m_m) {
        for (auto &element : element_lj_s) { // j1 = s, j2 = l
            if (element.second == std::numeric_limits<real_t>::max()) {

                element.second = pow(-1, k) * sqrt((2*element.first.j[0]+1)*(2*element.first.j[1]+1)) *
                        gsl_sf_coupling_6j(2*0.5, 2*element.first.j[0], 2*element.first.l[0], 2*element.first.j[1], 2*0.5, 2*k);
                //TODO: query
            }
        }
    }

    for (auto &element : element_lj_l) { // j1 = l, j2 = s
        if (element.second == std::numeric_limits<real_t>::max()) {

            element.second = pow(-1, k) * sqrt((2*element.first.j[0]+1)*(2*element.first.j[1]+1)) *
                    gsl_sf_coupling_6j(2*element.first.l[0], 2*element.first.j[0], 2*0.5, 2*element.first.j[1], 2*element.first.l[1], 2*k);
            //TODO: query
        }
    }

    for (auto &element : element_jm) {
        if (element.second == std::numeric_limits<real_t>::max()) {
            int q = element.first.m[0]-element.first.m[1];

            element.second = gsl_sf_coupling_3j(2*element.first.j[0], 2*k, 2*element.first.j[1], -2*element.first.m[0], 2*q, 2*element.first.m[1]);
            //TODO: query
        }
    }



    std::cout << 8.32 << std::endl;

    // write missing elements to database
    //TODO
}

real_t MatrixElements::getDipole(StateOne state_row, StateOne state_col) {
    return pow(-1, state_row.j-state_row.m+0.5+state_col.l+state_row.j) *
            element_nlj[StateTwo({{state_row.n, state_col.n}}, {{state_row.l, state_col.l}}, {{0,0}}, {{state_row.j, state_col.j}}, {{0,0}}).order()] *
            element_jm[StateTwo({{0,0}}, {{0, 0}}, {{0,0}}, {{state_row.j, state_col.j}}, {{state_row.m, state_col.m}}).order()] *
            element_lj_l[StateTwo({{0,0}}, {{state_row.l, state_col.l}}, {{0,0}}, {{state_row.j, state_col.j}}, {{0,0}}).order()];

}

real_t MatrixElements::getMomentum(StateOne state_row, StateOne state_col) {
    /*return pow(-1, state_row.j-state_row.m+state_row.l+0.5+state_col.j) * muB *
            element_jm[StateTwo({{0,0}}, {{0, 0}}, {{0,0}}, {{state_row.j, state_col.j}}, {{state_row.m, state_col.m}}).order()] *
            (gL*element_lj_l[StateTwo({{0,0}}, {{state_row.l, state_col.l}}, {{0,0}}, {{state_row.j, state_col.j}}, {{0,0}}).order()] *
            sqrt(state_row.l*(state_row.l+1)*(2*state_row.l+1)) +
            gS*element_lj_s[StateTwo({{0,0}}, {{state_row.l, state_col.l}}, {{0,0}}, {{state_row.j, state_col.j}}, {{0,0}}).order()] *
            sqrt(state_row.s*(state_row.s+1)*(2*state_row.s+1)));*/

    return pow(-1, state_row.j-state_row.m) * muB *
            element_jm[StateTwo({{0,0}}, {{0, 0}}, {{0,0}}, {{state_row.j, state_col.j}}, {{state_row.m, state_col.m}}).order()] *
            (gL*element_lj_l[StateTwo({{0,0}}, {{state_row.l, state_col.l}}, {{0,0}}, {{state_row.j, state_col.j}}, {{0,0}}).order()] *
            sqrt(state_row.l*(state_row.l+1)*(2*state_row.l+1)) * pow(-1, 0.5+state_col.l+state_row.j) +
            gS*element_lj_s[StateTwo({{0,0}}, {{state_row.l, state_col.l}}, {{0,0}}, {{state_row.j, state_col.j}}, {{0,0}}).order()] *
            sqrt(state_row.s*(state_row.s+1)*(2*state_row.s+1)) * pow(-1, 0.5+state_col.j+state_row.l));

    /*if ((state_row.j != state_col.j) || (state_row.l != state_col.l) || (state_row.m != state_col.m) || (state_row.n != state_col.n)) {
        real_t gJ = 3./2.+(state_row.s*(state_row.s+1)-state_row.l*(state_row.l+1))/(2*state_row.j*(state_row.j+1));
        return gJ*muB*state_row.m;

        return pow(-1, state_row.j-state_row.m+state_row.l+0.5+state_col.j) * muB *
                element_jm[StateTwo({{0,0}}, {{0, 0}}, {{0,0}}, {{state_row.j, state_col.j}}, {{state_row.m, state_col.m}}).order()] *
                (gS-gL)*element_lj_s[StateTwo({{0,0}}, {{state_row.l, state_col.l}}, {{0,0}}, {{state_row.j, state_col.j}}, {{0,0}}).order()] *
                sqrt(state_row.s*(state_row.s+1)*(2*state_row.s+1));
    }*/
}
