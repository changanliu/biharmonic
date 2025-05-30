using GeneralGraphs, GraphDatasets
using TOML, ProgressMeter, ProgressBars

function download_graph(d::AbstractDict)::NormalUnweightedGraph
    source, file_name = d["source"], d["file_name"]
    @assert source in ["SNAP", "KONECT"]
    @info "Details of graph $(d["name"])" source file_name n = d["total_n"] m = d["total_m"] LCC_n = d["LCC_n"] LCC_m = d["LCC_m"]
    if source == "SNAP"
        g = loadUndiSNAP(d["url"], d["name"])
        LCC!(g)
        NormalUnweightedGraph(g)
    else
        # loadUndiKONECT(d["internal_name"], d["name"]) |> LCC |> NormalUnweightedGraph
        g = loadUndiKONECT(d["internal_name"], d["name"])
        LCC!(g)
        NormalUnweightedGraph(g)
    end
end

function load_graph(d::AbstractDict, dir_path::AbstractString="data/input")::NormalUnweightedGraph
    source, file_name = d["source"], d["file_name"]
    file_path = joinpath(dir_path, "$file_name.txt")
    if !isfile(file_path)
        g = download_graph(d)
        write(file_path, g)
        g
    else
        @info "Details of graph $(d["name"])" source file_name n = d["total_n"] m = d["total_m"] LCC_n = d["LCC_n"] LCC_m = d["LCC_m"]
        NormalUnweightedGraph(d["name"], file_path)
    end
end

function download_data(d::AbstractDict, dir_path::AbstractString="data/input")
    source, file_name = d["source"], d["file_name"]
    file_path = joinpath(dir_path, "$file_name.txt")
    @assert source in ["SNAP", "KONECT"]
    if !isfile(file_path)
        write(file_path, download_graph(d))
    else
        @info "File of graph $(d["name"]) already exists"
    end
end