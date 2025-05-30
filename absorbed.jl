using Random
using LinearAlgebra
using Base.Threads

function write_list_to_csv(file::IO, list)
    str_list = join(string.(list), ",")
    write(file, str_list * "\n")
end

function node_pairs(n::Int, v::Int, len::Int)
    ans_set = Set{Tuple{Int,Int}}()
    nodes = Int[]
    sizehint!(nodes, n - 1)
    for u in 1:n
        if u != v
            push!(nodes, u)
        end
    end
    while length(ans_set) < len
        s, t = rand(nodes), rand(nodes)
        if s > t
            s, t = t, s
        end
        if s == t || (s, t) in ans_set
            continue
        end
        push!(ans_set, (s, t))
    end
    ans_set |> collect
end

function random_walk_until_absorb(path, l, nbr, start, absorb, lmax)
    l = l + 1
    current_node = start
    path[l] = current_node

    while current_node != absorb && l < lmax
        neighbors = nbr[current_node]
        current_node = neighbors[rand(1:end)]
        l = l + 1
        path[l] = current_node
    end

    return l
end

function calculate_variable(path1, path2, l1, l2, b1, b2, invdegsq, counts)
    value = 0.0
    for node in path2[b2:(l2 - 1)]
        counts[node] = counts[node] + 1
    end

    for node in path1[b1:(l1 - 1)]
        if counts[node] != 0
            invdegreesq = invdegsq[node]
            value += counts[node] * invdegreesq
        end
    end

    counts[path2[b2:l2]] .= 0
    return value
end
function special_random_walk(path, nbr, v, lmax)
    neighbors = nbr[v]
    l = 1
    first_step = neighbors[rand(1:end)]
    path[l] = v
    l = random_walk_until_absorb(path, l, nbr, first_step, v, lmax)
    return l
end

function simulate_random_walks(lmax, invdegsq, nbr, s, t, v, T)
    W, X, Y, Z, U, V, D = 0,0,0,0,0,0,0
    counts = zeros(Int, length(nbr))
    S1 = zeros(Int, lmax)
    S2 = zeros(Int, lmax)
    T1 = zeros(Int, lmax)
    T2 = zeros(Int, lmax)
    V1 = zeros(Int, lmax)
    V2 = zeros(Int, lmax)
    for i=1:T
        s1 = random_walk_until_absorb(S1, 0, nbr, s, v, lmax)
        s2 = random_walk_until_absorb(S2, 0, nbr, s, v, lmax)
        t1 = random_walk_until_absorb(T1, 0, nbr, t, v, lmax)
        t2 = random_walk_until_absorb(T2, 0, nbr, t, v, lmax)
        v1 = special_random_walk(V1, nbr, v, lmax)
        v2 = special_random_walk(V2, nbr, v, lmax)

        W += calculate_variable(S1, S2, s1, s2, 1, 1, invdegsq, counts)
        X += calculate_variable(T1, T2, t1, t2, 1, 1, invdegsq, counts)
        Y += calculate_variable(S1, T2, s1, t2, 1, 1, invdegsq, counts)
        Z += calculate_variable(S2, T1, s2, t1, 1, 1, invdegsq, counts)

        U += calculate_variable(S1, V1, s1, v1 - 1, 1, 2, invdegsq, counts)
        V += calculate_variable(T1, V2, t1, v2 - 1, 1, 2, invdegsq, counts)
        D += calculate_variable(V1, V2, v1 - 1, v2 - 1, 2, 2, invdegsq, counts)
    end
    W /= T; X /= T; Y /= T; Z /= T; U /= T; V /= T; D /= T
    return W, X, Y, Z, U, V, D
end

include("logw.jl")
include("graph.jl")
fname = open("./filename.txt", "r")
str = readline(fname);
nn = parse(Int, str);

filenames = []
for iii = 1:nn 
    str = readline(fname);
    str = strip(str);
    if '#' in str
        continue
    end
    push!(filenames, str)
end

@threads for iii = 1:length(filenames)
    str = filenames[iii]
    G = get_graph(str)
    n = G.n
    m = G.m
    L = lapsp(G);
    A = adjsp(G);

    dm = spzeros(n,n)
    deg = zeros(Int, n)
    invdegsq = zeros(n)
    for i=1:n
        dm[i,i] = L[i,i]
        deg[i] = L[i,i]
        invdegsq[i] = 1.0 / (deg[i] * deg[i])
    end
    v = argmax(deg)

    test_num = 20
    lmax = n
    test_pairs = node_pairs(n, v, test_num)

    @threads for T in [5000, 10000]
        errorfile = open(string(string("output/", str, "_$T",".csv")),"w")
        for (s, t) in test_pairs
            t1 = time()
            W, X, Y, Z, U, V, D = simulate_random_walks(lmax, invdegsq, G.nbr, s, t, v, T)
            APP = W + X - Y - Z - (deg[v] * (U - V))^2 / (1 + deg[v]^2 * D)
            t2 = time()
            write_list_to_csv(errorfile, [s, t, APP, t2 - t1])
        end
        close(errorfile)
    end
end