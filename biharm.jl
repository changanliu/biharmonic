using GeneralGraphs, GraphDatasets, LinearAlgebraUtils
using LinearAlgebra, SparseArrays, Base.Threads
using ProgressMeter, DataStructures

function spantree_biharm_parnew(g::NormalUnweightedGraph, v::Int, s::Int, t::Int, sample_num::Int)
    n = num_nodes(g)
    vis = falses(n)
    parent = Vector{Int}(undef, n)
    bfn = Int[]
    sizehint!(bfn, n - 1)
    q = Deque{Int}()
    push!(q, v)
    while !isempty(q)
        src = popfirst!(q)
        for u in g.adjs[src]
            if !vis[u]
                vis[u] = true
                push!(q, u)
                parent[u] = src
                push!(bfn, u)
            end
        end
    end

    t_num = nthreads()
    tot_in_forest = Matrix{Bool}(undef, n, t_num)
    tot_fa = Matrix{Int}(undef, n, t_num)
    tot_node_cur_b = Matrix{Int}(undef, n, t_num)
    tot_node_cur_r = Matrix{Int}(undef, n, t_num)
    tot_edge_cur_b = [Dict{Tuple{Int,Int},Int}() for _ in 1:t_num]
    tot_edge_cur_r = [Dict{Tuple{Int,Int},Int}() for _ in 1:t_num]
    @threads for t_id in 1:t_num
        for u in bfn
            f = parent[u]
            tot_edge_cur_b[t_id][(u, f)] = tot_edge_cur_r[t_id][(u, f)] = 0
        end
    end
    tot_chains = [Int[] for _ in 1:t_num]
    tot_chain = [Int[] for _ in 1:t_num]
    pm = Progress(sample_num; dt=0.5)
    @threads for _ in 1:sample_num
        t_id = threadid()
        in_forest = view(tot_in_forest, :, t_id)
        fa = view(tot_fa, :, t_id)
        node_cur_b = view(tot_node_cur_b, :, t_id)
        node_cur_r = view(tot_node_cur_r, :, t_id)
        edge_cur_b = tot_edge_cur_b[t_id]
        edge_cur_r = tot_edge_cur_r[t_id]
        chains = tot_chains[t_id]
        chain = tot_chain[t_id]

        fill!(in_forest, false)
        in_forest[v] = true
        empty!(chains)
        fill!(node_cur_b, 0)
        fill!(node_cur_r, 0)
        node_cur_b[s], node_cur_b[t] = 1, -1
        proc_nodes = Int[s, t]
        sizehint!(proc_nodes, length(g.adjs[v]) + 2)
        for u in g.adjs[v]
            push!(proc_nodes, u)
            node_cur_r[u] = -1
        end
        for src in proc_nodes
            u = src
            while !in_forest[u]
                fa[u] = rand(g.adjs[u])
                u = fa[u]
            end
            u = src
            empty!(chain)
            while !in_forest[u]
                in_forest[u] = true
                push!(chain, u)
                u = fa[u]
            end
            reverse!(chain)
            append!(chains, chain)
        end
        reverse!(chains)
        for u in chains
            f = fa[u]
            node_cur_b[f] += node_cur_b[u]
            node_cur_r[f] += node_cur_r[u]
            if parent[u] == f
                edge_cur_b[(u, f)] += node_cur_b[u]
                edge_cur_r[(u, f)] += node_cur_r[u]
            elseif parent[f] == u
                edge_cur_b[(f, u)] -= node_cur_b[u]
                edge_cur_r[(f, u)] -= node_cur_r[u]
            end
        end
        next!(pm)
    end
    finish!(pm)

    vecb, vecr = zeros(Int, n), zeros(Int, n)
    for u in bfn
        f = parent[u]
        vecb[u] = vecb[f] + sum(tot_edge_cur_b[t_id][(u, f)] for t_id in 1:t_num)
        vecr[u] = vecr[f] + sum(tot_edge_cur_r[t_id][(u, f)] for t_id in 1:t_num)
    end
    vecb = subvec(vecb / sample_num, v)
    vecr = subvec(vecr / sample_num, v)
    diagb = my_dot(vecb)
    denom = 1 + my_dot(vecr)
    numer = dot(vecb, vecr)^2
    diagb - numer / denom
end

function spantree_biharm_par(g::NormalUnweightedGraph, v::Int, s::Int, t::Int, sample_num::Int)
    n = num_nodes(g)
    vis = falses(n)
    parent = Vector{Int}(undef, n)
    bfn = Int[]
    sizehint!(bfn, n - 1)
    q = Deque{Int}()
    push!(q, v)
    while !isempty(q)
        src = popfirst!(q)
        for u in g.adjs[src]
            if !vis[u]
                vis[u] = true
                push!(q, u)
                parent[u] = src
                push!(bfn, u)
            end
        end
    end

    t_num = nthreads()
    tot_in_forest = Matrix{Bool}(undef, n, t_num)
    tot_fa = Matrix{Int}(undef, n, t_num)
    tot_node_cur_b = Matrix{Int}(undef, n, t_num)
    tot_node_cur_r = Matrix{Int}(undef, n, t_num)
    tot_edge_cur_b = [Dict{Tuple{Int,Int},Int}() for _ in 1:t_num]
    tot_edge_cur_r = [Dict{Tuple{Int,Int},Int}() for _ in 1:t_num]
    @threads for t_id in 1:t_num
        for u in bfn
            f = parent[u]
            tot_edge_cur_b[t_id][(u, f)] = tot_edge_cur_r[t_id][(u, f)] = 0
        end
    end
    tot_chains = [Int[] for _ in 1:t_num]
    tot_chain = [Int[] for _ in 1:t_num]
    pm = Progress(sample_num; dt=0.5)
    @threads for _ in 1:sample_num
        t_id = threadid()
        in_forest = view(tot_in_forest, :, t_id)
        fa = view(tot_fa, :, t_id)
        node_cur_b = view(tot_node_cur_b, :, t_id)
        node_cur_r = view(tot_node_cur_r, :, t_id)
        edge_cur_b = tot_edge_cur_b[t_id]
        edge_cur_r = tot_edge_cur_r[t_id]
        chains = tot_chains[t_id]
        chain = tot_chain[t_id]

        fill!(in_forest, false)
        in_forest[v] = true
        empty!(chains)
        fill!(node_cur_b, 0)
        fill!(node_cur_r, 0)
        node_cur_b[s], node_cur_b[t] = 1, -1
        foreach(u -> node_cur_r[u] = -1, g.adjs[v])
        for src in 1:n
            u = src
            while !in_forest[u]
                fa[u] = rand(g.adjs[u])
                u = fa[u]
            end
            u = src
            empty!(chain)
            while !in_forest[u]
                in_forest[u] = true
                push!(chain, u)
                u = fa[u]
            end
            reverse!(chain)
            append!(chains, chain)
        end
        reverse!(chains)
        for u in chains
            f = fa[u]
            node_cur_b[f] += node_cur_b[u]
            node_cur_r[f] += node_cur_r[u]
            if parent[u] == f
                edge_cur_b[(u, f)] += node_cur_b[u]
                edge_cur_r[(u, f)] += node_cur_r[u]
            elseif parent[f] == u
                edge_cur_b[(f, u)] -= node_cur_b[u]
                edge_cur_r[(f, u)] -= node_cur_r[u]
            end
        end
        next!(pm)
    end
    finish!(pm)

    vecb, vecr = zeros(Int, n), zeros(Int, n)
    for u in bfn
        f = parent[u]
        vecb[u] = vecb[f] + sum(tot_edge_cur_b[t_id][(u, f)] for t_id in 1:t_num)
        vecr[u] = vecr[f] + sum(tot_edge_cur_r[t_id][(u, f)] for t_id in 1:t_num)
    end
    vecb = subvec(vecb / sample_num, v)
    vecr = subvec(vecr / sample_num, v)
    diagb = my_dot(vecb)
    denom = 1 + my_dot(vecr)
    numer = dot(vecb, vecr)^2
    diagb - numer / denom
end
