struct wGraph
    n :: Int # |V|
    m :: Int # |E|
    nbr :: Array{Array{Tuple{Int, Float64}, 1}, 1} # neighbors of each vertex
end

struct Graph
    n :: Int # |V|
    m :: Int # |E|
    nbr :: Array{Array{Int, 1}, 1} # neighbors of each vertex
end

using Statistics

function rand_edge(G :: Graph)
    x = rand(1:G.n)
    y = G.nbr[x][rand(1:length(G.nbr[x]))]
    return (x, y)
end

# Graph -> Graph
function bfs_max_cc(G :: Graph)
    vis = fill(false, G.n)

    bfs(u :: Int) = begin
        vis[u] = true
        ret = [u]
        hinge = [u]
        while !isempty(hinge)
            v = pop!(hinge)
            for t in G.nbr[v]
                if !vis[t]
                    push!(ret, t)
                    push!(hinge,t)
                    vis[t] = true
                end
            end
        end
        return ret
    end
    max_cc = Vector{Int}(undef, 0)
    for u in 1 : G.n
        if ! vis[u]
            cc = bfs(u)
            if length(cc) > length(max_cc)
                max_cc = cc
            end
        end
    end

    id = Dict{Int, Int}()
    for i in eachindex(max_cc)
        id[max_cc[i]] = i
    end

    n = length(id)

    nbr = Array{Array{Int, 1}, 1}()
    for i in 1:n push!(nbr, []) end
    for u in filter(x -> haskey(id, x), eachindex(G.nbr))
        for v in filter(x -> haskey(id, x), G.nbr[u])
            push!(nbr[id[u]], id[v])
        end
    end

    m = div(sum(length.(nbr)), 2)

    return Graph(n, m, nbr)
end

# String -> Graph
function read_graph(str :: String)
    ints = parse.(Int, split(str))

    n = 0
    d = Dict{Int, Int}()
    edges = Set{Tuple{Int, Int}}()
    getid(x :: Int) = haskey(d, x) ? d[x] : d[x] = n += 1

    for i in 1 : 2 : length(ints)
        u = getid(ints[i])
        v = getid(ints[i + 1])
        if u == v continue end
        if u > v  u, v = v, u end
        push!(edges, (u, v))
    end
    m = length(edges)

    nbr = Array{Array{Int, 1}, 1}()
    for i in 1:n push!(nbr, []) end
    for (u, v) in edges
        push!(nbr[u], v)
        push!(nbr[v], u)
    end

    return Graph(n, m, nbr)
end

# String -> Graph
function read_file(filename :: String)
    return open(filename) do f
        #read_graph( read(f, String) )
        bfs_max_cc( read_graph( read(f, String) ) )
    end
end
