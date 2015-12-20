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
    species = "Rb";

    muB = 0.5;
    gS = 2.0023192;
    gL = 1;    
}

void MatrixElements::precalculate(BasisnamesOne basis_one, bool exist_d_0, bool exist_d_p, bool exist_d_m, bool exist_m_0, bool exist_m_p, bool exist_m_m) {
    //TODO: check if k = 1 in case of exist_m

    SQLite3 db("cache_matrix_elements.db");

    // create cache tables if necessary

    db.exec("CREATE TABLE IF NOT EXISTS cache_nlj ("
            "species text, k integer, n1 integer, l1 integer, j1 double,"
            "n2 integer, l2 integer, j2 double, value double, UNIQUE (species, k, n1, l1, j1, n2, l2, j2));");

    db.exec("CREATE TABLE IF NOT EXISTS cache_lj_s ("
            "species text, k integer, l1 integer, j1 double,"
            "l2 integer, j2 double, value double, UNIQUE (species, k, l1, j1, l2, j2));");


    db.exec("CREATE TABLE IF NOT EXISTS cache_lj_l ("
            "species text, k integer, l1 integer, j1 double,"
            "l2 integer, j2 double, value double, UNIQUE (species, k, l1, j1, l2, j2));");


    db.exec("CREATE TABLE IF NOT EXISTS cache_jm ("
            "species text, k integer, j1 integer, m1 double,"
            "j2 integer, m2 double, value double, UNIQUE (species, k, j1, m1, j2, m2));");

    // determine elements

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

                StateTwo state;

                if ((exist_d_0 || exist_d_p || exist_d_m) && (abs(state_row.l-state_col.l) == 1)) {
                    state = StateTwo({{state_row.n, state_col.n}}, {{state_row.l, state_col.l}}, {{0,0}}, {{state_row.j, state_col.j}}, {{0,0}}).order();
                    auto result_nlj = element_nlj.insert({state, std::numeric_limits<real_t>::max()});

                    if (result_nlj.second) {
                        ss.str(std::string());
                        ss << "insert into tmp_nlj (n1,l1,j1,n2,l2,j2) values ("
                           << state.n[0] << "," << state.l[0] << "," << state.j[0] << ","
                           << state.n[1] << "," << state.l[1] << "," << state.j[1] << ");";
                        db.exec(ss.str());
                    }
                }

                if ((exist_m_0 || exist_m_p || exist_m_m) && (state_row.l == state_col.l) && (state_row.n == state_col.n) ) {
                    state = StateTwo({{0,0}}, {{state_row.l, state_col.l}}, {{0,0}}, {{state_row.j, state_col.j}}, {{0,0}}).order();
                    auto result_lj = element_lj_s.insert({state, std::numeric_limits<real_t>::max()});

                    if (result_lj.second) {
                        ss.str(std::string());
                        ss << "insert into tmp_lj_s (l1,j1,l2,j2) values ("
                           << state.l[0] << "," << state.j[0] << ","
                           << state.l[1] << "," << state.j[1] << ");";
                        db.exec(ss.str());
                    }
                }

                state = StateTwo({{0,0}}, {{state_row.l, state_col.l}}, {{0,0}}, {{state_row.j, state_col.j}}, {{0,0}}).order();
                auto result_lj = element_lj_l.insert({state, std::numeric_limits<real_t>::max()});

                if (result_lj.second) {
                    ss.str(std::string());
                    ss << "insert into tmp_lj_l (l1,j1,l2,j2) values ("
                       << state.l[0] << "," << state.j[0] << ","
                       << state.l[1] << "," << state.j[1] << ");";
                    db.exec(ss.str());
                }

                state = StateTwo({{0,0}}, {{0, 0}}, {{0,0}}, {{state_row.j, state_col.j}}, {{state_row.m, state_col.m}}).order();
                auto result_jm = element_jm.insert({state, std::numeric_limits<real_t>::max()});

                if (result_jm.second) {
                    ss.str(std::string());
                    ss << "insert into tmp_jm (j1,m1,j2,m2) values ("
                       << state.j[0] << "," << state.m[0] << ","
                       << state.j[1] << "," << state.m[1] << ");";
                    db.exec(ss.str());
                }
            }
        }
    }

    db.exec("end transaction;");

    std::cout << 8.31 << std::endl;

    // load from database
    int n1, n2, l1, l2;
    float j1, j2, m1, m2;
    double value;

    if (exist_d_0 || exist_d_p || exist_d_m) {
        ss.str(std::string());
        ss << "SELECT c.n1, c.l1, c.j1, c.n2, c.l2, c.j2, c.value FROM cache_nlj c INNER JOIN tmp_nlj t ON ("
           << "c.n1 = t.n1 AND c.l1 = t.l1 AND c.j1 = t.j1 AND c.n2 = t.n2 AND c.l2 = t.l2 AND c.j2 = t.j2) "
           << "WHERE c.species = '" << species << "' AND c.k = " << k << ";";
        SQLite3Result result_nlj = db.query(ss.str());
        for (auto r : result_nlj) {
            r >> n1 >> l1 >> j1 >> n2 >> l2 >> j2 >> value;
            element_nlj[StateTwo({{n1, n2}}, {{l1, l2}}, {{0,0}}, {{j1, j2}}, {{0,0}})] = value;
        }
    }

    if (exist_m_0 || exist_m_p || exist_m_m) {
        ss.str(std::string());
        ss << "SELECT c.l1, c.j1, c.l2, c.j2, c.value FROM cache_lj_s c INNER JOIN tmp_lj_s t ON ("
           << "c.l1 = t.l1 AND c.j1 = t.j1 AND c.l2 = t.l2 AND c.j2 = t.j2) "
           << "WHERE c.species = '" << species << "' AND c.k = " << k << ";";
        SQLite3Result result_lj_s = db.query(ss.str());
        for (auto r : result_lj_s) {
            r >> l1 >> j1 >> l2 >> j2 >> value;
            element_lj_s[StateTwo({{0, 0}}, {{l1, l2}}, {{0,0}}, {{j1, j2}}, {{0,0}})] = value;
        }
    }

    ss.str(std::string());
    ss << "SELECT c.l1, c.j1, c.l2, c.j2, c.value FROM cache_lj_l c INNER JOIN tmp_lj_l t ON ("
       << "c.l1 = t.l1 AND c.j1 = t.j1 AND c.l2 = t.l2 AND c.j2 = t.j2) "
       << "WHERE c.species = '" << species << "' AND c.k = " << k << ";";
    SQLite3Result result_lj_l = db.query(ss.str());
    for (auto r : result_lj_l) {
        r >> l1 >> j1 >> l2 >> j2 >> value;
        element_lj_l[StateTwo({{0, 0}}, {{l1, l2}}, {{0,0}}, {{j1, j2}}, {{0,0}})] = value;
    }

    ss.str(std::string());
    ss << "SELECT c.j1, c.m1, c.j2, c.m2, c.value FROM cache_jm c INNER JOIN tmp_jm t ON ("
       << "c.j1 = t.j1 AND c.m1 = t.m1 AND c.j2 = t.j2 AND c.m2 = t.m2) "
       << "WHERE c.species = '" << species << "' AND c.k = " << k << ";";
    SQLite3Result result_jm = db.query(ss.str());
    for (auto r : result_jm) {
        r >> j1 >> m1 >> j2 >> m2 >> value;
        element_jm[StateTwo({{0, 0}}, {{0, 0}}, {{0,0}}, {{j1, j2}}, {{m1,m2}})] = value;
    }


    std::cout << 8.311 << std::endl;


    // calculate missing elements

    db.exec("begin transaction;");

    if (exist_d_0 || exist_d_p || exist_d_m) {
        for (auto &element : element_nlj) {
            if (element.second == std::numeric_limits<real_t>::max()) {
                int lmax = fmax(element.first.l[0],element.first.l[1]);

                element.second = pow(-1, element.first.l[0]+lmax) * sqrt(lmax) *
                        radial_element(species, element.first.n[0], element.first.l[0], element.first.j[0], k, element.first.n[1], element.first.l[1], element.first.j[1]); // TODO <Rb|Cs> not possible, so Rb is enough

                ss.str(std::string());
                ss << "insert into cache_nlj (species, k, n1, l1, j1, n2, l2, j2, value) values ("
                   << "'" << species << "'" << "," << k << ","
                   << element.first.n[0] << "," << element.first.l[0] << "," << element.first.j[0] << ","
                   << element.first.n[1] << "," << element.first.l[1] << "," << element.first.j[1] << ","
                   << element.second << ");";
                db.exec(ss.str());
            }
        }
    }

    std::cout << 8.312 << std::endl;

    if (exist_m_0 || exist_m_p || exist_m_m) {
        for (auto &element : element_lj_s) { // j1 = s, j2 = l
            if (element.second == std::numeric_limits<real_t>::max()) {

                element.second = pow(-1, k) * sqrt((2*element.first.j[0]+1)*(2*element.first.j[1]+1)) *
                        gsl_sf_coupling_6j(2*0.5, 2*element.first.j[0], 2*element.first.l[0], 2*element.first.j[1], 2*0.5, 2*k);

                ss.str(std::string());
                ss << "insert into cache_lj_s (species, k, l1, j1, l2, j2, value) values ("
                   << "'" << species << "'" << "," << k << ","
                   << element.first.l[0] << "," << element.first.j[0] << ","
                   << element.first.l[1] << "," << element.first.j[1] << ","
                   << element.second << ");";
                db.exec(ss.str());
            }
        }
    }

    for (auto &element : element_lj_l) { // j1 = l, j2 = s
        if (element.second == std::numeric_limits<real_t>::max()) {

            element.second = pow(-1, k) * sqrt((2*element.first.j[0]+1)*(2*element.first.j[1]+1)) *
                    gsl_sf_coupling_6j(2*element.first.l[0], 2*element.first.j[0], 2*0.5, 2*element.first.j[1], 2*element.first.l[1], 2*k);

            ss.str(std::string());
            ss << "insert into cache_lj_l (species, k, l1, j1, l2, j2, value) values ("
               << "'" << species << "'" << "," << k << ","
               << element.first.l[0] << "," << element.first.j[0] << ","
               << element.first.l[1] << "," << element.first.j[1] << ","
               << element.second << ");";
            db.exec(ss.str());
        }
    }

    for (auto &element : element_jm) {
        if (element.second == std::numeric_limits<real_t>::max()) {
            int q = element.first.m[0]-element.first.m[1];

            element.second = gsl_sf_coupling_3j(2*element.first.j[0], 2*k, 2*element.first.j[1], -2*element.first.m[0], 2*q, 2*element.first.m[1]);

            ss.str(std::string());
            ss << "insert into cache_jm (species, k, j1, m1, j2, m2, value) values ("
               << "'" << species << "'" << "," << k << ","
               << element.first.j[0] << "," << element.first.m[0] << ","
               << element.first.j[1] << "," << element.first.m[1] << ","
               << element.second << ");";
            db.exec(ss.str());
        }
    }

    db.exec("end transaction;");



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
