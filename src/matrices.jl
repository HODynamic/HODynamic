using LinearAlgebra
using SparseArrays

include("graph.jl")


function getAsp(G :: Graph)
	I, J, V :: Array{Float32, 1} = [], [], []
    for i in eachindex(G.nbr)
        for j in G.nbr[i]
            push!(I, i)
            push!(J, j)
            push!(V, 1)
        end
    end
    return sparse(I, J, V, G.n, G.n)
end


function getPsp(G :: wGraph)
	I, J, V :: Array{Float32, 1} = [], [], []
    for i in eachindex(G.nbr)
        for t in G.nbr[i]
            push!(I, i)
            push!(J, t[1])
            push!(V, t[2])
        end
    end
    A = sparse(I, J, V, G.n, G.n)
    dinv = Array{Float32, 1}(undef, G.n)
	for i in 1:G.n
		sum = 0.0
		for t in G.nbr[i]
			sum += t[2]
		end
		dinv[i] = 1.0 / sum
	end
    Dinv = spdiagm(0 => dinv)
    P = Dinv * A
    return P
end


function getP(G :: Graph)
	A = zeros(Float32, G.n, G.n)
	dinv = Array{Float32, 1}(undef, G.n)
	for i in 1:G.n
		dinv[i] = 1.0 / convert(Float32, size(G.nbr[i], 1))
	end
    Dinv = Diagonal(dinv)
    for i in eachindex(G.nbr)
        A[CartesianIndex.(tuple.(i, G.nbr[i]))] .=  1.0
    end
    P = Dinv * A
    return P
end


function getPast(G :: Graph, beta)
	P = getP(G :: Graph)
	curP = P
	Past = zeros(Float32, G.n, G.n)
	for b in beta
		# println(b, curP)
		Past += b * curP
		curP *= P
	end
	return Past
end

