#include "MatrixElements.h"
#include "SQLite.hpp"
#include <sstream>
#include <iostream>
#include <string>

bool selectionRulesDipole(StateOne state1, StateOne state2, int q) {
    return (abs(state1.l-state2.l) == 1) && (state1.m == state2.m+q) && (fabs(state1.j-state2.j) <= 1);
}

bool selectionRulesQuadrupole(StateOne state1, StateOne state2, int q) {
    return (abs(state1.l-state2.l) == 0 || abs(state1.l-state2.l) == 2) && (state1.m == state2.m+q) && (fabs(state1.j-state2.j) <= 2) && (state1.j != 0.5 || state2.j != 0.5);
}


bool selectionRulesMomentum(StateOne state1, StateOne state2, int q) {
    return (state1.l == state2.l) && (state1.m == state2.m+q) && (fabs(state1.j-state2.j) <= 1); // && (state1.n == state2.n);
}

size_t findidx(std::vector<real_t> x, real_t d) {
 size_t i;
 for (i = 0; i < x.size(); ++i) {
   if (x[i] == d)
     break;
 }
 return i;
}


MatrixElements::MatrixElements(std::string species, int k, std::string dbname) : species(species), k(k), dbname(dbname), calc_missing(true) {
    muB = 0.5;
    gS = 2.0023192;
    gL = 1;
}


MatrixElements::MatrixElements(const Configuration& config, std::string species, int k, std::string dbname) : MatrixElements(species,k,dbname) {
    calc_missing = config["calcmissing"].str() == "true";
}

void MatrixElements::precalculate_momentum(std::shared_ptr<const BasisnamesOne> basis_one, bool exist_0, bool exist_p, bool exist_m) {
    precalculate(basis_one,false,false,false,false,false,false,false,false,exist_0,exist_p,exist_m);
}

void MatrixElements::precalculate_dipole(std::shared_ptr<const BasisnamesOne> basis_one, bool exist_0, bool exist_p, bool exist_m) {
    precalculate(basis_one,exist_0,exist_p,exist_m,false,false,false,false,false,false,false,false);
}

void MatrixElements::precalculate_quadrupole(std::shared_ptr<const BasisnamesOne> basis_one, bool exist_0, bool exist_p, bool exist_m, bool exist_pp, bool exist_mm) {
    precalculate(basis_one,false,false,false,exist_0,exist_p,exist_m,exist_pp,exist_mm,false,false,false);
}

void MatrixElements::precalculate(std::shared_ptr<const BasisnamesOne> basis_one, bool exist_d_0, bool exist_d_p, bool exist_d_m, bool exist_q_0, bool exist_q_p, bool exist_q_m, bool exist_q_pp, bool exist_q_mm, bool exist_m_0, bool exist_m_p, bool exist_m_m) {
    if (!(exist_m_0 || exist_m_p || exist_m_m || exist_d_0 || exist_d_p || exist_d_m || exist_q_0 || exist_q_p || exist_q_m || exist_q_pp || exist_q_mm)) return;

    if ((exist_m_0 || exist_m_p || exist_m_m) && k != 1) {
        std::cout << "For calculating momentum matrix elements, it must be k=1." << std::endl;
        abort();
    }

    if ((exist_d_0 || exist_d_p || exist_d_m) && k != 1) {
        std::cout << "For calculating dipole matrix elements, it must be k=1." << std::endl;
        abort();
    }

    if ((exist_q_0 || exist_q_p || exist_q_m || exist_q_pp || exist_q_mm) && k != 2) {
        std::cout << "For calculating quadrupole matrix elements, it must be k=2." << std::endl;
        abort();
    }

    SQLite3 db(dbname);

    // create cache tables if necessary

    db.exec("CREATE TABLE IF NOT EXISTS cache_nlj_k ("
            "species text, k integer, n1 integer, l1 integer, j1 double,"
            "n2 integer, l2 integer, j2 double, value double, UNIQUE (species, k, n1, l1, j1, n2, l2, j2));");

    db.exec("CREATE TABLE IF NOT EXISTS cache_nlj_0 ("
            "species text, n1 integer, l1 integer, j1 double,"
            "n2 integer, l2 integer, j2 double, value double, UNIQUE (species, n1, l1, j1, n2, l2, j2));");

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

    db.exec("CREATE TEMPORARY TABLE tmp_nlj_k ("
            "n1 integer, l1 integer, j1 double,"
            "n2 integer, l2 integer, j2 double);");

    db.exec("CREATE TEMPORARY TABLE tmp_nlj_0 ("
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

    for (const auto &state_col : *basis_one) {
        for (const auto &state_row : *basis_one) {
            if (state_row.idx < state_col.idx) {
                continue;
            }

            if ((exist_d_0 && selectionRulesDipole(state_row, state_col, 0)) ||
                    (exist_d_p && selectionRulesDipole(state_row, state_col, 1)) ||
                    (exist_d_m && selectionRulesDipole(state_row, state_col, -1)) ||
                    (exist_q_0 && selectionRulesQuadrupole(state_row, state_col, 0)) ||
                    (exist_q_p && selectionRulesQuadrupole(state_row, state_col, 1)) ||
                    (exist_q_m && selectionRulesQuadrupole(state_row, state_col, -1)) ||
                    (exist_q_pp && selectionRulesQuadrupole(state_row, state_col, 2)) ||
                    (exist_q_mm && selectionRulesQuadrupole(state_row, state_col, -2)) ||
                    (exist_m_0 && selectionRulesMomentum(state_row, state_col, 0)) ||
                    (exist_m_p && selectionRulesMomentum(state_row, state_col, 1)) ||
                    (exist_m_m && selectionRulesMomentum(state_row, state_col, -1))) {

                StateTwo state;

                if (((exist_d_0 || exist_d_p || exist_d_m) && (abs(state_row.l-state_col.l) == 1)) ||
                        ((exist_q_0 || exist_q_p || exist_q_m || exist_q_pp || exist_q_mm) && (abs(state_row.l-state_col.l) == 0 || abs(state_row.l-state_col.l) == 2))) {
                    state = StateTwo({{state_row.n, state_col.n}}, {{state_row.l, state_col.l}}, {{0,0}}, {{state_row.j, state_col.j}}, {{0,0}}).order();
                    auto result_nlj_k = element_nlj_k.insert({state, std::numeric_limits<real_t>::max()});

                    if (result_nlj_k.second) {
                        ss.str(std::string());
                        ss << "insert into tmp_nlj_k (n1,l1,j1,n2,l2,j2) values ("
                           << state.n[0] << "," << state.l[0] << "," << state.j[0] << ","
                           << state.n[1] << "," << state.l[1] << "," << state.j[1] << ");";
                        db.exec(ss.str());
                    }
                }

                if ((exist_m_0 || exist_m_p || exist_m_m) && (state_row.l == state_col.l)) {
                    state = StateTwo({{0,0}}, {{state_row.l, state_col.l}}, {{0,0}}, {{state_row.j, state_col.j}}, {{0,0}}).order();
                    auto result_lj = element_lj_s.insert({state, std::numeric_limits<real_t>::max()});

                    if (result_lj.second) {
                        ss.str(std::string());
                        ss << "insert into tmp_lj_s (l1,j1,l2,j2) values ("
                           << state.l[0] << "," << state.j[0] << ","
                           << state.l[1] << "," << state.j[1] << ");";
                        db.exec(ss.str());
                    }

                    state = StateTwo({{state_row.n, state_col.n}}, {{state_row.l, state_col.l}}, {{0,0}}, {{state_row.j, state_col.j}}, {{0,0}}).order();
                    auto result_nlj_0 = element_nlj_0.insert({state, std::numeric_limits<real_t>::max()});

                    if (result_nlj_0.second) {
                        ss.str(std::string());
                        ss << "insert into tmp_nlj_0 (n1,l1,j1,n2,l2,j2) values ("
                           << state.n[0] << "," << state.l[0] << "," << state.j[0] << ","
                           << state.n[1] << "," << state.l[1] << "," << state.j[1] << ");";
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

    // load from database
    int n1, n2, l1, l2;
    float j1, j2, m1, m2;
    double value;

    if ((exist_d_0 || exist_d_p || exist_d_m) || (exist_q_0 || exist_q_p || exist_q_m || exist_q_pp || exist_q_mm)) {
        ss.str(std::string());
        ss << "SELECT c.n1, c.l1, c.j1, c.n2, c.l2, c.j2, c.value FROM cache_nlj_k c INNER JOIN tmp_nlj_k t ON ("
           << "c.n1 = t.n1 AND c.l1 = t.l1 AND c.j1 = t.j1 AND c.n2 = t.n2 AND c.l2 = t.l2 AND c.j2 = t.j2) "
           << "WHERE c.species = '" << species << "' AND c.k = " << k << ";";
        SQLite3Result result_nlj_k = db.query(ss.str());
        for (auto r : result_nlj_k) {
            *r >> n1 >> l1 >> j1 >> n2 >> l2 >> j2 >> value;
            element_nlj_k[StateTwo({{n1, n2}}, {{l1, l2}}, {{0,0}}, {{j1, j2}}, {{0,0}})] = value;
        }
    }

    if (exist_m_0 || exist_m_p || exist_m_m) {
        ss.str(std::string());
        ss << "SELECT c.l1, c.j1, c.l2, c.j2, c.value FROM cache_lj_s c INNER JOIN tmp_lj_s t ON ("
           << "c.l1 = t.l1 AND c.j1 = t.j1 AND c.l2 = t.l2 AND c.j2 = t.j2) "
           << "WHERE c.species = '" << species << "' AND c.k = " << k << ";";
        SQLite3Result result_lj_s = db.query(ss.str());
        for (auto r : result_lj_s) {
            *r >> l1 >> j1 >> l2 >> j2 >> value;
            element_lj_s[StateTwo({{0, 0}}, {{l1, l2}}, {{0,0}}, {{j1, j2}}, {{0,0}})] = value;
        }

        ss.str(std::string());
        ss << "SELECT c.n1, c.l1, c.j1, c.n2, c.l2, c.j2, c.value FROM cache_nlj_0 c INNER JOIN tmp_nlj_0 t ON ("
           << "c.n1 = t.n1 AND c.l1 = t.l1 AND c.j1 = t.j1 AND c.n2 = t.n2 AND c.l2 = t.l2 AND c.j2 = t.j2) "
           << "WHERE c.species = '" << species << "';";
        SQLite3Result result_nlj_0 = db.query(ss.str());
        for (auto r : result_nlj_0) {
            *r >> n1 >> l1 >> j1 >> n2 >> l2 >> j2 >> value;
            element_nlj_0[StateTwo({{n1, n2}}, {{l1, l2}}, {{0,0}}, {{j1, j2}}, {{0,0}})] = value;
        }
    }

    ss.str(std::string());
    ss << "SELECT c.l1, c.j1, c.l2, c.j2, c.value FROM cache_lj_l c INNER JOIN tmp_lj_l t ON ("
       << "c.l1 = t.l1 AND c.j1 = t.j1 AND c.l2 = t.l2 AND c.j2 = t.j2) "
       << "WHERE c.species = '" << species << "' AND c.k = " << k << ";";
    SQLite3Result result_lj_l = db.query(ss.str());
    for (auto r : result_lj_l) {
        *r >> l1 >> j1 >> l2 >> j2 >> value;
        element_lj_l[StateTwo({{0, 0}}, {{l1, l2}}, {{0,0}}, {{j1, j2}}, {{0,0}})] = value;
    }

    ss.str(std::string());
    ss << "SELECT c.j1, c.m1, c.j2, c.m2, c.value FROM cache_jm c INNER JOIN tmp_jm t ON ("
       << "c.j1 = t.j1 AND c.m1 = t.m1 AND c.j2 = t.j2 AND c.m2 = t.m2) "
       << "WHERE c.species = '" << species << "' AND c.k = " << k << ";";
    SQLite3Result result_jm = db.query(ss.str());
    for (auto r : result_jm) {
        *r >> j1 >> m1 >> j2 >> m2 >> value;
        element_jm[StateTwo({{0, 0}}, {{0, 0}}, {{0,0}}, {{j1, j2}}, {{m1,m2}})] = value;
    }

    /*std::cout << "Same n" << std::endl;
    std::cout << calcRadialElement(species, 50, 0, 0.5, 0, 50, 0, 0.5)<< std::endl;
    std::cout << calcRadialElement(species, 50, 0, 0.5, 0, 50, 0, 1.5)<< std::endl;
    std::cout << calcRadialElement(species, 50, 0, 0.5, 0, 50, 1, 0.5)<< std::endl;

    std::cout << "Different n" << std::endl;
    std::cout << calcRadialElement(species, 50, 0, 0.5, 0, 51, 0, 0.5)<< std::endl;
    std::cout << calcRadialElement(species, 50, 0, 0.5, 0, 51, 0, 1.5)<< std::endl;
    std::cout << calcRadialElement(species, 50, 0, 0.5, 0, 51, 1, 0.5)<< std::endl;*/

    // calculate missing elements and write them to the database

    if (calc_missing) {

        db.exec("begin transaction;");

        if ((exist_d_0 || exist_d_p || exist_d_m) || (exist_q_0 || exist_q_p || exist_q_m || exist_q_pp || exist_q_mm)) {
            for (auto &element : element_nlj_k) {
                if (element.second == std::numeric_limits<real_t>::max()) {
                    int lmax = fmax(element.first.l[0],element.first.l[1]);

                    if (k == 1) { // dipole
                        element.second = pow(-1, lmax) * sqrt(lmax) *
                            calcRadialElement(species, element.first.n[0], element.first.l[0], element.first.j[0], k, element.first.n[1], element.first.l[1], element.first.j[1]);
                    } else if (k == 2) { // quadrupole
                        //element.second = 0; // TODO !!!
                        element.second = pow(-1, lmax) * sqrt(lmax) *
                            calcRadialElement(species, element.first.n[0], element.first.l[0], element.first.j[0], k, element.first.n[1], element.first.l[1], element.first.j[1]);
                    }

                    ss.str(std::string());
                    ss << "insert into cache_nlj_k (species, k, n1, l1, j1, n2, l2, j2, value) values ("
                       << "'" << species << "'" << "," << k << ","
                       << element.first.n[0] << "," << element.first.l[0] << "," << element.first.j[0] << ","
                       << element.first.n[1] << "," << element.first.l[1] << "," << element.first.j[1] << ","
                       << element.second << ");";
                    db.exec(ss.str());
                }
            }
        }

        if (exist_m_0 || exist_m_p || exist_m_m) {
            for (auto &element : element_lj_s) { // j1 = s, j2 = l
                if (element.second == std::numeric_limits<real_t>::max()) {

                    //element.second = pow(-1, k) * sqrt((2*element.first.j[0]+1)*(2*element.first.j[1]+1)) *
                    //        gsl_sf_coupling_6j(2*0.5, 2*element.first.j[0], 2*element.first.l[0], 2*element.first.j[1], 2*0.5, 2*k);

                    element.second = pow(-1, k) * sqrt((2*element.first.j[0]+1)*(2*element.first.j[1]+1)) *
                        WignerSymbols::wigner6j(0.5, element.first.j[0], element.first.l[0], element.first.j[1], 0.5, k);

                    ss.str(std::string());
                    ss << "insert into cache_lj_s (species, k, l1, j1, l2, j2, value) values ("
                       << "'" << species << "'" << "," << k << ","
                       << element.first.l[0] << "," << element.first.j[0] << ","
                       << element.first.l[1] << "," << element.first.j[1] << ","
                       << element.second << ");";
                    db.exec(ss.str());
                }
            }

            for (auto &element : element_nlj_0) {
                if (element.second == std::numeric_limits<real_t>::max()) {

                    element.second =
                        calcRadialElement(species, element.first.n[0], element.first.l[0], element.first.j[0], 0, element.first.n[1], element.first.l[1], element.first.j[1]);

                    ss.str(std::string());
                    ss << "insert into cache_nlj_0 (species, n1, l1, j1, n2, l2, j2, value) values ("
                       << "'" << species << "'" << ","
                       << element.first.n[0] << "," << element.first.l[0] << "," << element.first.j[0] << ","
                       << element.first.n[1] << "," << element.first.l[1] << "," << element.first.j[1] << ","
                       << element.second << ");";
                    db.exec(ss.str());
                }
            }
        }

        for (auto &element : element_lj_l) { // j1 = l, j2 = s
            if (element.second == std::numeric_limits<real_t>::max()) {

                //element.second = pow(-1, k) * sqrt((2*element.first.j[0]+1)*(2*element.first.j[1]+1)) *
                //        gsl_sf_coupling_6j(2*element.first.l[0], 2*element.first.j[0], 2*0.5, 2*element.first.j[1], 2*element.first.l[1], 2*k);

                element.second = pow(-1, k) * sqrt((2*element.first.j[0]+1)*(2*element.first.j[1]+1)) *
                    WignerSymbols::wigner6j(element.first.l[0], element.first.j[0], 0.5, element.first.j[1], element.first.l[1], k);

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

                //element.second = gsl_sf_coupling_3j(2*element.first.j[0], 2*k, 2*element.first.j[1], -2*element.first.m[0], 2*q, 2*element.first.m[1]);
                element.second = WignerSymbols::wigner3j(element.first.j[0], k, element.first.j[1], -element.first.m[0],q,element.first.m[1]);

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
    }
}

real_t MatrixElements::calcRadialElement(std::string species, int n1, int l1, real_t j1, int power,
                      int n2, int l2, real_t j2) {
  Numerov N1(species, n1, l1, j1);
  Numerov N2(species, n2, l2, j2);

  std::vector<real_t> x1 = N1.axis();
  std::vector<real_t> y1 = N1.integrate();
  std::vector<real_t> x2 = N2.axis();
  std::vector<real_t> y2 = N2.integrate();

  real_t xmin = N1.xmin >= N2.xmin ? N1.xmin : N2.xmin;
  real_t xmax = N1.xmax <= N2.xmax ? N1.xmax : N2.xmax;

  real_t mu = 0;
  // If there is an overlap, calculate the matrix element
  if (xmin <= xmax) {
    int start1 = findidx(x1, xmin);
    int end1   = findidx(x1, xmax);
    int start2 = findidx(x2, xmin);
    int end2   = findidx(x2, xmax);

    int i1, i2;
    for (i1 = start1, i2 = start2; i1 < end1 && i2 < end2; i1++, i2++) {
      mu += y1[i1]*y2[i2] * pow(x1[i1], 2*power+2) * N1.dx;
    }
    mu = fabs(2*mu);
  }

  return mu;
}

real_t MatrixElements::getDipole(StateOne state_row, StateOne state_col) {
    if (k != 1) abort();

    return pow(-1, state_row.j-state_row.m+0.5+state_col.l+state_row.j) * pow(-1, state_row.l) *
            element_nlj_k[StateTwo({{state_row.n, state_col.n}}, {{state_row.l, state_col.l}}, {{0,0}}, {{state_row.j, state_col.j}}, {{0,0}}).order()] *
            element_jm[StateTwo({{0,0}}, {{0, 0}}, {{0,0}}, {{state_row.j, state_col.j}}, {{state_row.m, state_col.m}}).order()] *
            element_lj_l[StateTwo({{0,0}}, {{state_row.l, state_col.l}}, {{0,0}}, {{state_row.j, state_col.j}}, {{0,0}}).order()];
}

real_t MatrixElements::getQuadrupole(StateOne state_row, StateOne state_col) {
    if (k != 2) abort();

    return pow(-1, state_row.j-state_row.m+0.5+state_col.l+state_row.j) * pow(-1, state_row.l) *
            element_nlj_k[StateTwo({{state_row.n, state_col.n}}, {{state_row.l, state_col.l}}, {{0,0}}, {{state_row.j, state_col.j}}, {{0,0}}).order()] * // TODO !!!
            element_jm[StateTwo({{0,0}}, {{0, 0}}, {{0,0}}, {{state_row.j, state_col.j}}, {{state_row.m, state_col.m}}).order()] *
            element_lj_l[StateTwo({{0,0}}, {{state_row.l, state_col.l}}, {{0,0}}, {{state_row.j, state_col.j}}, {{0,0}}).order()];
}

real_t MatrixElements::getMomentum(StateOne state_row, StateOne state_col) {
    return pow(-1, state_row.j-state_row.m) * muB *
            element_nlj_0[StateTwo({{state_row.n, state_col.n}}, {{state_row.l, state_col.l}}, {{0,0}}, {{state_row.j, state_col.j}}, {{0,0}}).order()] *
            element_jm[StateTwo({{0,0}}, {{0, 0}}, {{0,0}}, {{state_row.j, state_col.j}}, {{state_row.m, state_col.m}}).order()] *
            (gL*element_lj_l[StateTwo({{0,0}}, {{state_row.l, state_col.l}}, {{0,0}}, {{state_row.j, state_col.j}}, {{0,0}}).order()] *
            sqrt(state_row.l*(state_row.l+1)*(2*state_row.l+1)) * pow(-1, 0.5+state_col.l+state_row.j) +
            gS*element_lj_s[StateTwo({{0,0}}, {{state_row.l, state_col.l}}, {{0,0}}, {{state_row.j, state_col.j}}, {{0,0}}).order()] *
            sqrt(state_row.s*(state_row.s+1)*(2*state_row.s+1)) * pow(-1, 0.5+state_col.j+state_row.l));
}
