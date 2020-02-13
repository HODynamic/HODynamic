include("graph.jl")
include("matrices.jl")


function iterSolve(G :: wGraph, ALPHA, s, iters)
    P = getPsp(G)
    M = (I - ALPHA) * P
    z = copy(s)
    as = ALPHA * s
	for i in 1:iters
		z = as + M * z
	end
    return z
end
