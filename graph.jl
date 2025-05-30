struct Graph
    n :: Int # |V|
    m :: Int # |E|
    u :: Array{Int, 1}
    v :: Array{Int, 1} # uv is an edge
    w :: Array{Int, 1} # edge weight
#    w :: Array{Float64, 1} # weight of each edge
    nbr :: Array{Array{Int, 1}, 1}
    # Origin :: Dict{Int32, Int32} #origin label
    name :: String
    degree :: Array{Int, 1}
end


include("core.jl")

function get_model_graph(ffname)
    n = 0
    Label = Dict{Int32, Int32}()
    Origin = Dict{Int32, Int32}()

    getID(x :: Int) = haskey(Label, x) ? Label[x] : Label[x] = n += 1
    fname = string("./graphs/model/",ffname,".txt")
    fin = open(fname, "r")

    str = readline(fin)
    u = Int[]
    v = Int[]
    str = strip(str)
    tot = 0
    i=0
    while str != ""
        if occursin("#", str)
            str = strip(readline(fin))
            continue
        end
        str = split(str)
        if length(str)!=2
            str = strip(readline(fin))
            continue
        end
        i = i + 1
        # println(str)
        if length(str)<2
            println("error", i)
        end
        
        # println(str)
        x   = parse(Int, str[1])
        y   = parse(Int, str[2])
        # if x<y
        # println(x, " ", y)
        u1 = getID(x)
        v1 = getID(y)
        Origin[u1] = x
        Origin[v1] = y
        push!(u, u1)
        push!(v, v1)
        # println(tot)
        tot += 1
        # end
        str = strip(readline(fin))
    end
    nbr=[ [ ] for i in 1:n ]
    for i=1:tot
        u1=u[i];
        v1=v[i];
        push!(nbr[u1],v1);
        push!(nbr[v1],u1);
    end

    degree = zeros(n)
    for i=1:n
        degree[i] = length(nbr[i])
    end
    w = zeros(Float64, tot)
    for i=1:tot
        w[i] = 1
    end

    close(fin)
    return Graph(n, tot, u, v, w, nbr, ffname, degree)
end

function get_graph(ffname)
    n = 0
    Label = Dict{Int32, Int32}()
    Origin = Dict{Int32, Int32}()

    getID(x :: Int) = haskey(Label, x) ? Label[x] : Label[x] = n += 1
    fname = string("./graphs/",ffname,".txt")
    fin = open(fname, "r")

    str = readline(fin)
    u = Int[]
    v = Int[]
    str = strip(str)
    tot = 0
    i=0
    while str != ""
        if occursin("#", str)
            str = strip(readline(fin))
            continue
        end
        str = split(str)
        if length(str)!=2
            str = strip(readline(fin))
            continue
        end
        i = i + 1
        # println(str)
        if length(str)<2
            println("error", i)
        end
        
        # println(str)
        x   = parse(Int, str[1])
        y   = parse(Int, str[2])
        # if x<y
        # println(x, " ", y)
        u1 = getID(x)
        v1 = getID(y)
        Origin[u1] = x
        Origin[v1] = y
        push!(u, u1)
        push!(v, v1)
        # println(tot)
        tot += 1
        # end
        str = strip(readline(fin))
    end
    nbr=[ [ ] for i in 1:n ]
    for i=1:tot
        u1=u[i];
        v1=v[i];
        push!(nbr[u1],v1);
        push!(nbr[v1],u1);
    end

    degree = zeros(n)
    for i=1:n
        degree[i] = length(nbr[i])
    end
    w = zeros(Float64, tot)
    for i=1:tot
        w[i] = 1
    end

    close(fin)
    return Graph(n, tot, u, v, w, nbr, ffname, degree)
end


function reorder_graph_by_degree(g::Graph)
    n = g.n

    idx = collect(1:n)
    
    sort_idx = sortperm(g.degree)
    
    old_to_new = Dict{Int, Int}()
    for i in 1:n
        old_to_new[sort_idx[i]] = i
    end
    
    new_u = Array{Int,1}(undef, g.m)
    new_v = Array{Int,1}(undef, g.m)
    new_nbr = [Int[] for i in 1:n]
    
    for i in 1:g.m
        new_u[i] = old_to_new[g.u[i]]
        new_v[i] = old_to_new[g.v[i]]
        push!(new_nbr[new_u[i]], new_v[i])
        push!(new_nbr[new_v[i]], new_u[i])
    end
    
    new_degree = zeros(Int, n)
    for i in 1:n
        new_degree[i] = length(new_nbr[i])
    end

    return Graph(n, g.m, new_u, new_v, g.w, new_nbr, g.name, new_degree)
end
