# julia main.jl exact/sparsify/approx/all Iter/M/beta/best uni/ex/pl

include("graph.jl")
include("logw.jl")
include("matrices.jl")
include("sparsify.jl")
include("IterSolve.jl")


function ex(n, xmin=1, Lambda=1)
	s = zeros(Float32, n, 1)
	for i in 1:n
		s[i] = xmin - (1.0/Lambda) * log(1-rand())
	end
	s ./= maximum(s)
	return s
end


function pl(n, alpha, xmin=1)
	s = zeros(Float32, n, 1)
	for i in 1:n
		s[i] = xmin * (1.0-rand())^(-1.0/(alpha-1.0))
	end
	s ./= maximum(s)
	return s
end


run_exact = (ARGS[1] == "exact" || ARGS[1] == "sparsify" || ARGS[1] == "all" )
run_sparsify = (ARGS[1] == "sparsify" || ARGS[1] == "all")
run_iter = (ARGS[1] == "approx" || ARGS[1] == "all")

datadir = string("../data/")
outFName = string("../", ARGS[1], "_", ARGS[2], "_", ARGS[3], ".txt")
w = open(outFName, "w")

(ARGS[2] == "beta") ? (betas = [[1.0f0], [0.0f0, 1.0f0]]) : (betas = [[0.5f0, 0.5f0]])
(ARGS[2] == "Iter") ? (canIter = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 50, 100, 200, 500, 1000]) : (canIter = [100])
(ARGS[2] == "M") ? (canM = [1, 10, 100, 200, 500, 1000, 2000]) : (canM = [10])

for rFile in filter(x -> endswith(x, ".txt"), readdir(string(datadir)))
    logw(w, "reading graph from edges list file ", rFile)
    G = read_file(string(datadir, rFile))
    logw(w, "finished reading graph")
    logw(w, "LCC.n: ", G.n, "\t LCC.m: ", G.m)
    flush(stdout)
	
	if ARGS[3] == "uni"
		s = rand(Float32, (G.n, 1))
	elseif ARGS[3] == "pl"
		s = pl(G.n, 2.5)
	else
		s = ex(G.n)
	end
	alpha = rand(Float32, G.n)
	# s = [1.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0]
	# alpha = [1.0; 0.6; 0.6; 0.6; 0.01; 0.01; 0.01; 0.01; 0.01; 0.01]
	ALPHA = Diagonal(alpha)
	#logw(w, "s: ", s)
	logw(w, "sum: ", sum(s))
	#logw(w, "alpha: ", alpha[1:10])
	zs = []
	for beta in betas
		logw(w, "computing solution for beta = ", beta)
		# Exact Solution
		if run_exact
			logw(w, "")
			logw(w, "computing exact solution")
			start_time = time()
			Past = getPast(G, beta)
			logw(w, "get Past")
			z_1 = inv(I - (I - ALPHA) * Past) * ALPHA * s
			elapsed_time = time() - start_time
			logw(w, "finished computation")
			logw(w, "z: ", z_1[1:3])
			logw(w, "sum: ", sum(z_1))
			logw(w, "time: ", elapsed_time)
			push!(zs, z_1)
		end

		# Sparisify
		for M in canM
			logw(w, "")
	    	logw(w, "start sparsification for M = ", M)
	    	start_time = time()
	    	Gtil = getGtil(G, beta, M)
	    	elapsed_time = time() - start_time
			logw(w, "finished sparsification")
			logw(w, "Gtil.n: ", Gtil.n, "\t Gtil.m: ", Gtil.m)
			logw(w, "time: ", elapsed_time)
		
			# Exact Solution on Sparsifier
			if run_sparsify
				logw(w, "")
	    		logw(w, "computing exact solution on sparsifier")
	    		start_time = time()
				Ptil = Array(getPsp(Gtil))
				z_2 = inv(I - (I - ALPHA) * Ptil) * ALPHA * s
				elapsed_time = time() - start_time
				logw(w, "finished computation")
				logw(w, "z: ", z_2[1:3])
				logw(w, "time: ", elapsed_time)
				if run_exact 
					logw(w, "max error: ", maximum(abs.(z_1-z_2)))
					logw(w, "avg error: ", sum(abs.(z_1-z_2)/G.n)) 
				end
			end
	
			# Approximate Solution via Iteration
			if run_iter
				for iter in canIter
					logw(w, "")
	    			logw(w, "computing approximate solution via Iteration for Iter = ", iter)
	    			start_time = time()
	    			z_3 = iterSolve(Gtil, ALPHA, s, iter)
					elapsed_time = time() - start_time
					logw(w, "finished computation")
					logw(w, "z: ", z_3[1:3])
					logw(w, "time: ", elapsed_time)
					if run_exact
						logw(w, "max error: ", maximum(abs.(z_1-z_3)))
						logw(w, "avg error: ", sum(abs.(z_1-z_3)/G.n))
					end
				end
			end
		end
	end
	
	for i in 2:size(zs, 1)
		sp = [0f0, 0.01f0, 0.05f0, 0.1f0, 0.15f0, 0.2f0, 1f0]
		for j in 2:size(sp,1)
			l, r = sp[j-1], sp[j]
			function cv(x)
				return (r>=x>l) ? 1 : 0
			end
			logw(w, "difference in ($l, $r] : ", sum(cv.(abs.(zs[i]-zs[1])))/G.n * 100, '%')
		end
	end
	
    logw(w, "")
    logw(w, String(fill('*', 80)))
    logw(w, "")
end

close(w)
