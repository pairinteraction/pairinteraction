# TestQuantumDefect
using Test
using SQLite
using PairInteraction

# test_comparison
qd = PairInteraction.QuantumDefect("Rb", 78, 1, 0.5)

sqldb = SQLite.DB("pairinteraction/databases/quantum_defects.db")
res = SQLite.Query(sqldb, 
                   "select ac,Z,a1,a2,a3,a4,rc from model_potential where ( (element = 'Rb') and (L = 1) );")

res_tuple, = [res...]

@test PairInteraction.ac(qd) == res_tuple.ac
@test PairInteraction.Z(qd) == res_tuple.Z
@test PairInteraction.a1(qd) == res_tuple.a1
@test PairInteraction.a2(qd) == res_tuple.a2
@test PairInteraction.a3(qd) == res_tuple.a3
@test PairInteraction.a4(qd) == res_tuple.a4
@test PairInteraction.rc(qd) == res_tuple.rc
