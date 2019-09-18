module PairInteraction

import Base.^
import Base.string
using CxxWrap

pairinteraction_dir = dirname(@__FILE__)
build_dir = joinpath(splitpath(pairinteraction_dir)[1:end-2]..., "build", "pairinteraction")

@wrapmodule(joinpath(build_dir, "libpireal.so"))

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
