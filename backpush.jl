include("logw.jl")
include("graph.jl")
include("core.jl")
using Base.Threads
using DataStructures
using TOML

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


function vbackpush(G::Graph, v::Int, s::Int, t::Int, eps::Float64, r_max::Float64)
    push_threshold = r_max
    n = G.n
    deg = G.degree
    nbr = G.nbr

    r_s = zeros(Float64, n)
    r_t = zeros(Float64, n)
    r_v = zeros(Float64, n)

    r_s[s] = 1
    r_t[t] = 1

    for neighbor in nbr[v]
        r_v[neighbor] = 1 / deg[neighbor]
    end

    tau_v_s = zeros(Float64, n)
    tau_v_t = zeros(Float64, n)
    tau_v = zeros(Float64, n)

    r_sum = 1
    isInQueue = falses(n)
    r_q = Queue{Int}()
    enqueue!(r_q, s)
    isInQueue[s] = true
    while !isempty(r_q)
        current_node = dequeue!(r_q)
        isInQueue[current_node] = false

        tau_v_s[current_node] += r_s[current_node]
        increment = r_s[current_node]
        r_s[current_node] = 0.0

        for neighbor in G.nbr[current_node]
            if neighbor == v
                r_sum -= increment / deg[current_node]
            else
                update_node = neighbor
                r_s[update_node] += increment / deg[update_node]

                if !isInQueue[update_node] && r_s[update_node] >= (push_threshold * deg[s])
                    enqueue!(r_q, update_node)
                    isInQueue[update_node] = true
                end
            end
        end
    end

    r_sum = 1
    isInQueue = falses(n)
    r_q = Queue{Int}()
    enqueue!(r_q, t)
    isInQueue[t] = true
    while !isempty(r_q)
        current_node = dequeue!(r_q)
        isInQueue[current_node] = false
        tau_v_t[current_node] += r_t[current_node]
        increment = r_t[current_node]
        r_t[current_node] = 0.0
        for neighbor in G.nbr[current_node]
            if neighbor == v
                r_sum -= increment / deg[current_node]
            else
                update_node = neighbor
                r_t[update_node] += increment / deg[update_node]
                if !isInQueue[update_node] && r_t[update_node] >= (push_threshold * deg[t])
                    enqueue!(r_q, update_node)
                    isInQueue[update_node] = true
                end
            end
        end
    end

    isInQueue = falses(n)
    r_q = Queue{Int}()
    for neighbor in nbr[v]
        enqueue!(r_q, neighbor)
        isInQueue[neighbor] = true
    end
    while !isempty(r_q)
        current_node = dequeue!(r_q)
        isInQueue[current_node] = false
        tau_v[current_node] += r_v[current_node]
        increment = r_v[current_node]
        r_v[current_node] = 0.0

        for neighbor in G.nbr[current_node]
            if neighbor != v
                update_node = neighbor
                r_v[update_node] += increment / deg[update_node]
                
                if !isInQueue[update_node] && r_v[update_node] >= (push_threshold)
                    enqueue!(r_q, update_node)
                    isInQueue[update_node] = true
                end
            end
        end
    end

    X, Y, Z = 0.0, 0.0, 0.0
    for u in 1:n
        if u != v
        X += (tau_v_s[u]/deg[s] - tau_v_t[u]/deg[t])^2
        Y += (tau_v_s[u]/deg[s] - tau_v_t[u]/deg[t]) * tau_v[u]
        Z += tau_v[u]^2
        end
    end

    beta_prime = X - (Y^2) / (1 + Z)
    return beta_prime, X, Y^2, 1 + Z
end

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
    G0 = get_graph(str)
    G = reorder_graph_by_degree(G0)
    n = G.n
    m = G.m
    L = lapsp(G);

    dm = spzeros(n,n)
    deg = zeros(n)
    for i=1:n
        dm[i,i] = L[i,i]
        deg[i] = L[i,i]
    end
    st_idx = sortperm(deg, rev=true);
    v = st_idx[1]
    @show v
    test_num = 100
    test_pairs = node_pairs(n, v, test_num)
    for thres in [0.0001, 0.00001, 0.000001]
        errorfile = open(string("output/", str, "_backpush_$thres",".csv"), "w")
        for (s, t) in test_pairs
            tt1 = time()
            biharmonic_distance, q11, q22, q33 = vbackpush(G, v, s, t, 0.0, thres)
            tt2 = time()
            write_list_to_csv(errorfile, [s, t, biharmonic_distance, tt2 - tt1])
        end
        close(errorfile)
    end
end
