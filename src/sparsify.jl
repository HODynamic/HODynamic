using LinearAlgebra
using SparseArrays
# using StatsBase

include("graph.jl")
include("matrices.jl")


function walk(G :: Graph, u, k)
	for i in 1:k
		u = G.nbr[u][rand(1:size(G.nbr[u], 1))]
	end
	return u
end


function sparsify(G :: Graph, edgelist, M, t)
	mp = Dict()
	@time for i in 1:M
		e = edgelist[rand(1:size(edgelist, 1))]
		r = rand(1:t)
		u = walk(G, e[1], r-1)
		v = walk(G, e[2], t-r)
		Z = 2 * t
		w = (2 * t * G.m) / (M * Z)
		haskey(mp, (u, v)) ? mp[(u, v)] += w : mp[(u, v)] = w
		(u != v) && (haskey(mp, (v, u)) ? mp[(v, u)] += w : mp[(v, u)] = w)
	end
	Is, Js, Vs = Array{Int, 1}(), Array{Int, 1}(), Array{Float32, 1}()
	@time for k in keys(mp)
		push!(Is, k[1])
		push!(Js, k[2])
		push!(Vs, mp[k])
	end
	@time A = sparse(Is, Js, Vs, G.n, G.n)
	return A
end


function getGtil(G :: Graph, beta, k)
	A = spzeros(Float32, G.n, G.n)
	edgelist = []
	for i in eachindex(G.nbr)
        for t in G.nbr[i]
            push!(edgelist, (i, t[1]))
        end
    end
	for i in eachindex(beta)
		if i == 1
			Gi = getAsp(G)
		else
			M = k * i * G.m
			Gi = sparsify(G, edgelist, M, i)
		end
		A = A + beta[i] * Gi
	end
	
	n = G.n
	m = 0
	nbr = Array{Array{Tuple{Int, Float64}, 1}, 1}()
    for i in 1:n push!(nbr, []) end
    u, v, w = findnz(A)
	for i in eachindex(u)
		push!(nbr[u[i]], (v[i], w[i])) 
		(v[i] >= u[i]) && (m += 1)
	end
	return wGraph(n, m, nbr)
end

