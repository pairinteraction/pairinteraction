# Copyright (c) 2020 Sebastian Weber, Henri Menke, Alexander Papageorge. All rights reserved.
# 
# This file is part of the pairinteraction library.
#
# The pairinteraction library is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# The pairinteraction library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with the pairinteraction library. If not, see <http://www.gnu.org/licenses/>.

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
