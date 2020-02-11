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

module PairInteraction

import Base.^
import Base.string
using CxxWrap

pairinteraction_dir = dirname(@__FILE__)
@wrapmodule(joinpath(pairinteraction_dir, "libpireal_jl.so"))

function __init__()
    @initcxx
end

# Extend sparse.
import SparseArrays.sparse
function sparse(x::PairInteraction.eigen_sparse_tRef)
    rows = 1 .+ copy(PairInteraction.innerIndex(x))
    outer_starts = 1 .+ copy(PairInteraction.outerIndex(x))

    push!(outer_starts, length(rows))
    cols = zeros(Int, length(rows))

    for jdx in 1:length(outer_starts)-1
        cols[outer_starts[jdx]:outer_starts[jdx+1]] .= jdx
    end

    values = PairInteraction.nonzerorealvalues(x) + 1im * PairInteraction.nonzeroimagvalues(x)
    return sparse(rows, cols, values)
end
export sparse

end  # End module.
