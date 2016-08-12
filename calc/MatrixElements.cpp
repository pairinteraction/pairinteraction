#include "MatrixElements.h"
#include "SQLite.hpp"
#include <sstream>
#include <iostream>
#include <string>
#include <limits>

bool selectionRulesMomentum(StateOne state1, StateOne state2, int q) {
    bool validL = state1.l == state2.l;
    bool validJ = fabs(state1.j-state2.j) <= 1;
    bool validM = state1.m == state2.m+q;
    bool validQ = abs(q) <= 1;
    return validL && validJ && validM && validQ;
}

bool selectionRulesMomentum(StateOne state1, StateOne state2) {
    bool validL = state1.l == state2.l;
    bool validJ = fabs(state1.j-state2.j) <= 1;
    bool validM = (fabs(state1.m - state2.m) <= 1);
    return validL && validJ && validM;
}

bool selectionRulesMultipole(StateOne state1, StateOne state2, int kappa, int q) {
    bool validL = (abs(state1.l-state2.l) <= kappa) && (kappa%2 == abs(state1.l-state2.l)%2);
    bool validJ = (fabs(state1.j-state2.j) <= kappa)  && (state1.j+state2.j >= kappa);
    bool validM = state1.m == state2.m+q;
    bool validQ = abs(q) <= kappa;
    bool noZero = !(kappa == 2 && state1.j == state2.j && state2.j == 1.5 && state1.m == -state2.m && fabs(state1.m - state2.m) == 1);
    return validL && validJ && validM && validQ && noZero;
}

bool selectionRulesMultipole(StateOne state1, StateOne state2, int kappa) {
    bool validL = (abs(state1.l-state2.l) <= kappa) && (kappa%2 == abs(state1.l-state2.l)%2);
    bool validJ = (fabs(state1.j-state2.j) <= kappa)  && (state1.j+state2.j >= kappa);
    bool validM = (fabs(state1.m - state2.m) <= kappa);
    bool noZero = !(kappa == 2 && state1.j == state2.j && state2.j == 1.5 && state1.m == -state2.m && fabs(state1.m - state2.m) == 1);
    return validL && validJ && validM && noZero;
}


MatrixElements::MatrixElements(std::string species, std::string dbname) : species(species), dbname(dbname) {
    muB = 0.5;
    gS = 2.0023192;
    gL = 1;
}

MatrixElements::MatrixElements(const Configuration& config, std::string species, std::string dbname) : MatrixElements(species,dbname) {
    if (config["missingCalc"].str() == "true") {
        method = "Modelpotentials";
    } else if (config["missingWhittaker"].str() == "true") {
        method = "Whittaker";
    } else {
        method = "Error";
    }
}

void MatrixElements::precalculateMultipole(std::shared_ptr<const BasisnamesOne> basis_one, int k) {
    int m = std::numeric_limits<int>::max();
    precalculate(basis_one, k, m, k, true, false);
}

void MatrixElements::precalculateElectricMomentum(std::shared_ptr<const BasisnamesOne> basis_one, int q) {
    precalculate(basis_one, 1, q, 1, true, false);
}

void MatrixElements::precalculateMagneticMomentum(std::shared_ptr<const BasisnamesOne> basis_one, int q) {
    precalculate(basis_one, 1, q, 0, false, true);
}

void MatrixElements::precalculateDiamagnetism(std::shared_ptr<const BasisnamesOne> basis_one, int k, int q) {
    precalculate(basis_one, k, q, 2, true, false);
}

real_t MatrixElements::getElectricMomentum(StateOne state_row, StateOne state_col) {
    return getMultipole(state_row,state_col,1);
}

real_t MatrixElements::getMagneticMomentum(StateOne state_row, StateOne state_col) {
    real_t val =  muB * cache_radial[0][StateTwo({{state_row.n, state_col.n}}, {{state_row.l, state_col.l}}, {{0,0}}, {{state_row.j, state_col.j}}, {{0,0}}).order()] *
            cache_angular[1][StateTwo({{0, 0}}, {{0, 0}}, {{0,0}}, {{state_row.j, state_col.j}}, {{state_row.m, state_col.m}})] *
            (gL * cache_reduced_commutes_s[1][StateTwo({{0, 0}}, {{state_row.l, state_col.l}}, {{0,0}}, {{state_row.j, state_col.j}}, {{0, 0}})] *
            sqrt(state_row.l*(state_row.l+1)*(2*state_row.l+1))  +
            gS * cache_reduced_commutes_l[1][StateTwo({{0, 0}}, {{state_row.l, state_col.l}}, {{0,0}}, {{state_row.j, state_col.j}}, {{0, 0}})] *
            sqrt(state_row.s*(state_row.s+1)*(2*state_row.s+1)));
    return val;
}

real_t MatrixElements::getDiamagnetism(StateOne state_row, StateOne state_col, int k) {
    real_t val =    1./12.*cache_radial[2][StateTwo({{state_row.n, state_col.n}}, {{state_row.l, state_col.l}}, {{0,0}}, {{state_row.j, state_col.j}}, {{0,0}}).order()] *
            cache_angular[k][StateTwo({{0, 0}}, {{0, 0}}, {{0,0}}, {{state_row.j, state_col.j}}, {{state_row.m, state_col.m}})] *
            cache_reduced_commutes_s[k][StateTwo({{0, 0}}, {{state_row.l, state_col.l}}, {{0,0}}, {{state_row.j, state_col.j}}, {{0, 0}})] *
            cache_reduced_multipole[k][StateTwo({{0, 0}}, {{state_row.l, state_col.l}}, {{0,0}}, {{0, 0}}, {{0, 0}})];
    return val;
}

real_t MatrixElements::getMultipole(StateOne state_row, StateOne state_col, int k) {
    real_t val =  cache_radial[k][StateTwo({{state_row.n, state_col.n}}, {{state_row.l, state_col.l}}, {{0,0}}, {{state_row.j, state_col.j}}, {{0,0}}).order()] *
            cache_angular[k][StateTwo({{0, 0}}, {{0, 0}}, {{0,0}}, {{state_row.j, state_col.j}}, {{state_row.m, state_col.m}})] *
            cache_reduced_commutes_s[k][StateTwo({{0, 0}}, {{state_row.l, state_col.l}}, {{0,0}}, {{state_row.j, state_col.j}}, {{0, 0}})] *
            cache_reduced_multipole[k][StateTwo({{0, 0}}, {{state_row.l, state_col.l}}, {{0,0}}, {{0, 0}}, {{0, 0}})];
    return val;
}

void MatrixElements::precalculate(std::shared_ptr<const BasisnamesOne> basis_one, int kappa, int q, int kappar, bool calcElectricMultipole, bool calcMagneticMomentum) {
    SQLite3 db(dbname);

    // --- create cache tables if necessary (reduced_moemntumS and reduced_moemntumL need not to be cached since they are trivial) ---

    if (calcElectricMultipole || calcMagneticMomentum) {
        db.exec("CREATE TABLE IF NOT EXISTS cache_radial ("
                "method text, species text, k integer, n1 integer, l1 integer, j1 double,"
                "n2 integer, l2 integer, j2 double, value double, UNIQUE (method, species, k, n1, l1, j1, n2, l2, j2));");
        db.exec("CREATE TEMPORARY TABLE tmp_radial ("
                "n1 integer, l1 integer, j1 double,"
                "n2 integer, l2 integer, j2 double);");
    }

    if (calcElectricMultipole || calcMagneticMomentum) {
        db.exec("CREATE TABLE IF NOT EXISTS cache_angular ("
                "k integer, j1 double, m1 double,"
                "j2 double, m2 double, value double, UNIQUE (k, j1, m1, j2, m2));");
        db.exec("CREATE TEMPORARY TABLE tmp_angular ("
                "j1 double, m1 double,"
                "j2 double, m2 double);");
    }

    if (calcElectricMultipole || calcMagneticMomentum) {
        db.exec("CREATE TABLE IF NOT EXISTS cache_reduced_commutes_s ("
                "k integer, l1 integer, j1 double,"
                "l2 integer, j2 double, value double, UNIQUE (k, l1, j1, l2, j2));");
        db.exec("CREATE TEMPORARY TABLE tmp_reduced_commutes_s ("
                "l1 integer, j1 double,"
                "l2 integer, j2 double);");
    }

    if (calcMagneticMomentum) {
        db.exec("CREATE TABLE IF NOT EXISTS cache_reduced_commutes_l ("
                "k integer, l1 integer, j1 double,"
                "l2 integer, j2 double, value double, UNIQUE (k, l1, j1, l2, j2));");
        db.exec("CREATE TEMPORARY TABLE tmp_reduced_commutes_l ("
                "l1 integer, j1 double,"
                "l2 integer, j2 double);");
    }

    if (calcElectricMultipole) {
        db.exec("CREATE TABLE IF NOT EXISTS cache_reduced_multipole ("
                "k integer, l1 integer,"
                "l2 integer, value double, UNIQUE (k, l1, l2));");
        db.exec("CREATE TEMPORARY TABLE tmp_reduced_multipole ("
                "l1 integer,"
                "l2 integer);");
    }

    // --- determine elements ---

    std::stringstream ss;

    db.exec("begin transaction;");

    for (const auto &state_col : *basis_one) {
        for (const auto &state_row : *basis_one) {
            if (q != std::numeric_limits<int>::max() && state_row.m-state_col.m != q) {
                continue;
            }

            //if (state_row.idx < state_col.idx) { // TODO
            //    continue;
            //}

            if ((calcElectricMultipole && selectionRulesMultipole(state_row, state_col, kappa)) || (calcMagneticMomentum && selectionRulesMomentum(state_row, state_col))) {
                if (calcElectricMultipole || calcMagneticMomentum) {
                    StateTwo state_nlj = StateTwo({{state_row.n, state_col.n}}, {{state_row.l, state_col.l}}, {{0,0}}, {{state_row.j, state_col.j}}, {{0,0}}).order();
                    auto missing_cache_radial = cache_radial[kappar].insert({state_nlj, std::numeric_limits<real_t>::max()});
                    if (missing_cache_radial.second) {
                        ss.str(std::string());
                        ss << "insert into tmp_radial (n1,l1,j1,n2,l2,j2) values ("
                           << state_nlj.n[0] << "," << state_nlj.l[0] << "," << state_nlj.j[0] << ","
                           << state_nlj.n[1] << "," << state_nlj.l[1] << "," << state_nlj.j[1] << ");";
                        db.exec(ss.str());
                    }
                }

                if (calcElectricMultipole || calcMagneticMomentum) {
                    StateTwo state_jm = StateTwo({{0, 0}}, {{0, 0}}, {{0,0}}, {{state_row.j, state_col.j}}, {{state_row.m, state_col.m}});
                    auto missing_cache_angular = cache_angular[kappa].insert({state_jm, std::numeric_limits<real_t>::max()});
                    if (missing_cache_angular.second) {
                        ss.str(std::string());
                        ss << "insert into tmp_angular (j1,m1,j2,m2) values ("
                           << state_jm.j[0] << "," << state_jm.m[0] << ","
                           << state_jm.j[1] << "," << state_jm.m[1] << ");";
                        db.exec(ss.str());
                    }
                }

                if (calcElectricMultipole || calcMagneticMomentum) {
                    StateTwo state_lj = StateTwo({{0, 0}}, {{state_row.l, state_col.l}}, {{0,0}}, {{state_row.j, state_col.j}}, {{0, 0}});
                    auto missing_cache_reduced_commutes_s = cache_reduced_commutes_s[kappa].insert({state_lj, std::numeric_limits<real_t>::max()});
                    if (missing_cache_reduced_commutes_s.second) {
                        ss.str(std::string());
                        ss << "insert into tmp_reduced_commutes_s (l1,j1,l2,j2) values ("
                           << state_lj.l[0] << "," << state_lj.j[0] << ","
                           << state_lj.l[1] << "," << state_lj.j[1] << ");";
                        db.exec(ss.str());
                    }
                }

                if (calcMagneticMomentum) {
                    StateTwo state_lj = StateTwo({{0, 0}}, {{state_row.l, state_col.l}}, {{0,0}}, {{state_row.j, state_col.j}}, {{0, 0}});
                    auto missing_cache_reduced_commutes_l = cache_reduced_commutes_l[kappa].insert({state_lj, std::numeric_limits<real_t>::max()});
                    if (missing_cache_reduced_commutes_l.second) {
                        ss.str(std::string());
                        ss << "insert into tmp_reduced_commutes_l (l1,j1,l2,j2) values ("
                           << state_lj.l[0] << "," << state_lj.j[0] << ","
                           << state_lj.l[1] << "," << state_lj.j[1] << ");";
                        db.exec(ss.str());
                    }
                }

                if (calcElectricMultipole) {
                    StateTwo state_l = StateTwo({{0, 0}}, {{state_row.l, state_col.l}}, {{0,0}}, {{0, 0}}, {{0, 0}});
                    auto missing_cache_reduced_multipole = cache_reduced_multipole[kappa].insert({state_l, std::numeric_limits<real_t>::max()});
                    if (missing_cache_reduced_multipole.second) {
                        ss.str(std::string());
                        ss << "insert into tmp_reduced_multipole (l1,l2) values ("
                           << state_l.l[0] << ","
                           << state_l.l[1] << ");";
                        db.exec(ss.str());
                    }
                }
            }
        }
    }

    db.exec("end transaction;");

    // --- load from database ---

    int n1, n2, l1, l2;
    float j1, j2, m1, m2;
    double value;

    if (calcElectricMultipole || calcMagneticMomentum) {
        ss.str(std::string());
        ss << "SELECT c.n1, c.l1, c.j1, c.n2, c.l2, c.j2, c.value FROM cache_radial c INNER JOIN tmp_radial t ON ("
           << "c.n1 = t.n1 AND c.l1 = t.l1 AND c.j1 = t.j1 AND c.n2 = t.n2 AND c.l2 = t.l2 AND c.j2 = t.j2) "
           << "WHERE c.method= '" << method << "' AND c.species = '" << species << "' AND c.k = " << kappar << ";";
        SQLite3Result result_radial = db.query(ss.str());
        for (auto r : result_radial) {
            *r >> n1 >> l1 >> j1 >> n2 >> l2 >> j2 >> value;
            cache_radial[kappar][StateTwo({{n1, n2}}, {{l1, l2}}, {{0,0}}, {{j1, j2}}, {{0,0}})] = value;
        }
    }

    if (calcElectricMultipole || calcMagneticMomentum) {
        ss.str(std::string());
        ss << "SELECT c.j1, c.m1, c.j2, c.m2, c.value FROM cache_angular c INNER JOIN tmp_angular t ON ("
           << "c.j1 = t.j1 AND c.m1 = t.m1 AND c.j2 = t.j2 AND c.m2 = t.m2) "
           << "WHERE c.k = " << kappa << ";";
        SQLite3Result result_angular = db.query(ss.str());
        for (auto r : result_angular) {
            *r >> j1 >> m1 >> j2 >> m2 >> value;
            cache_angular[kappa][StateTwo({{0, 0}}, {{0, 0}}, {{0,0}}, {{j1, j2}}, {{m1,m2}})] = value;
        }
    }

    if (calcElectricMultipole || calcMagneticMomentum) {
        ss.str(std::string());
        ss << "SELECT c.l1, c.j1, c.l2, c.j2, c.value FROM cache_reduced_commutes_s c INNER JOIN tmp_reduced_commutes_s t ON ("
           << "c.l1 = t.l1 AND c.j1 = t.j1 AND c.l2 = t.l2 AND c.j2 = t.j2) "
           << "WHERE c.k = " << kappa << ";";
        SQLite3Result result_reduced_commutes_s = db.query(ss.str());
        for (auto r : result_reduced_commutes_s) {
            *r >> l1 >> j1 >> l2 >> j2 >> value;
            cache_reduced_commutes_s[kappa][StateTwo({{0, 0}}, {{l1, l2}}, {{0,0}}, {{j1, j2}}, {{0,0}})] = value;
        }
    }

    if (calcMagneticMomentum) {
        ss.str(std::string());
        ss << "SELECT c.l1, c.j1, c.l2, c.j2, c.value FROM cache_reduced_commutes_l c INNER JOIN tmp_reduced_commutes_l t ON ("
           << "c.l1 = t.l1 AND c.j1 = t.j1 AND c.l2 = t.l2 AND c.j2 = t.j2) "
           << "WHERE c.k = " << kappa << ";";
        SQLite3Result result_reduced_commutes_l = db.query(ss.str());
        for (auto r : result_reduced_commutes_l) {
            *r >> l1 >> j1 >> l2 >> j2 >> value;
            cache_reduced_commutes_l[kappa][StateTwo({{0, 0}}, {{l1, l2}}, {{0,0}}, {{j1, j2}}, {{0,0}})] = value;
        }
    }

    if (calcElectricMultipole) {
        ss.str(std::string());
        ss << "SELECT c.l1, c.l2, c.value FROM cache_reduced_multipole c INNER JOIN tmp_reduced_multipole t ON ("
           << "c.l1 = t.l1 AND c.l2 = t.l2) "
           << "WHERE c.k = " << kappa << ";";
        SQLite3Result result_reduced_multipole = db.query(ss.str());
        for (auto r : result_reduced_multipole) {
            *r >> l1 >> l2 >> value;
            cache_reduced_multipole[kappa][StateTwo({{0, 0}}, {{l1, l2}}, {{0,0}}, {{0, 0}}, {{0,0}})] = value;
        }
    }

    // --- calculate missing elements and write them to the database ---

    db.exec("begin transaction;");

    if (calcElectricMultipole || calcMagneticMomentum) {
        for (auto &cache : cache_radial[kappar]) {
            if (cache.second == std::numeric_limits<real_t>::max()) {
                auto state = cache.first;

                cache.second = calcRadialElement(species, state.n[0], state.l[0], state.j[0], kappar, state.n[1], state.l[1], state.j[1]);

                ss.str(std::string());
                ss << "insert into cache_radial (method, species, k, n1, l1, j1, n2, l2, j2, value) values ("
                   << "'" << method << "'" << "," << "'" << species << "'" << "," << kappar << ","
                   << state.n[0] << "," << state.l[0] << "," << state.j[0] << ","
                   << state.n[1] << "," << state.l[1] << "," << state.j[1] << ","
                   << cache.second << ");";
                db.exec(ss.str());

            }
        }
    }

    if (calcElectricMultipole || calcMagneticMomentum) {
        for (auto &cache : cache_angular[kappa]) {
            if (cache.second == std::numeric_limits<real_t>::max()) {
                auto state = cache.first;
                float q = state.m[0]-state.m[1];

                cache.second = pow(-1, state.j[0]-state.m[0]) *
                        WignerSymbols::wigner3j(state.j[0], kappa, state.j[1], -state.m[0], q, state.m[1]);

                ss.str(std::string());
                ss << "insert into cache_angular (k, j1, m1, j2, m2, value) values ("
                   << kappa << ","
                   << state.j[0] << "," << state.m[0] << ","
                   << state.j[1] << "," << state.m[1] << ","
                   << cache.second << ");";
                db.exec(ss.str());
            }
        }
    }

    if (calcElectricMultipole || calcMagneticMomentum) {
        for (auto &cache : cache_reduced_commutes_s[kappa]) {
            if (cache.second == std::numeric_limits<real_t>::max()) {
                auto state = cache.first;

                cache.second = pow(-1, state.l[0] + 0.5 + state.j[1] + kappa) * sqrt((2*state.j[0]+1)*(2*state.j[1]+1)) *
                        WignerSymbols::wigner6j(state.l[0], state.j[0], 0.5, state.j[1], state.l[1], kappa);

                ss.str(std::string());
                ss << "insert into cache_reduced_commutes_s (k, l1, j1, l2, j2, value) values ("
                   << kappa << ","
                   << state.l[0] << "," << state.j[0] << ","
                   << state.l[1] << "," << state.j[1] << ","
                   << cache.second << ");";
                db.exec(ss.str());
            }
        }
    }

    if (calcMagneticMomentum) {
        for (auto &cache : cache_reduced_commutes_l[kappa]) {
            if (cache.second == std::numeric_limits<real_t>::max()) {
                auto state = cache.first;

                cache.second = pow(-1, state.l[0] + 0.5 + state.j[0] + kappa) * sqrt((2*state.j[0]+1)*(2*state.j[1]+1)) *
                        WignerSymbols::wigner6j(0.5, state.j[0], state.l[0], state.j[1], 0.5, kappa);

                ss.str(std::string());
                ss << "insert into cache_reduced_commutes_l (k, l1, j1, l2, j2, value) values ("
                   << kappa << ","
                   << state.l[0] << "," << state.j[0] << ","
                   << state.l[1] << "," << state.j[1] << ","
                   << cache.second << ");";
                db.exec(ss.str());
            }
        }
    }

    if (calcElectricMultipole) {
        for (auto &cache : cache_reduced_multipole[kappa]) {
            if (cache.second == std::numeric_limits<real_t>::max()) {
                auto state = cache.first;

                cache.second = pow(-1, state.l[0]) * sqrt((2*state.l[0]+1)*(2*state.l[1]+1)) *
                        WignerSymbols::wigner3j(state.l[0], kappa, state.l[1], 0, 0, 0);

                ss.str(std::string());
                ss << "insert into cache_reduced_multipole (k, l1, l2, value) values ("
                   << kappa << ","
                   << state.l[0] << ","
                   << state.l[1] << ","
                   << cache.second << ");";
                db.exec(ss.str());
            }
        }
    }

    db.exec("end transaction;");
}

real_t MatrixElements::calcRadialElement(std::string species, int n1, int l1, real_t j1, int power,
                                         int n2, int l2, real_t j2) {
    if (method == "Modelpotentials") {
        Numerov N1(species, n1, l1, j1);
        Numerov N2(species, n2, l2, j2);
        return IntegrateRadialElement(N1, power, N2);
    } else if(method == "Whittaker") {
        Whittaker N1(species, n1, l1, j1);
        Whittaker N2(species, n2, l2, j2);
        return IntegrateRadialElement(N1, power, N2);
    } else {
        std::cout << ">>ERR" << "You have to provide all radial matrix elements on your own because you have deactivated the calculation of missing radial matrix elements!" << std::endl; // TODO make it thread save
        abort();
    }
}
