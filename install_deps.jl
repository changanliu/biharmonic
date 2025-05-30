using Pkg

const urls = [
    "https://gitee.com/xhs7700/GeneralGraphs.jl#main",
    "https://gitee.com/xhs7700/GraphDatasets.jl#main",
    "https://gitee.com/xhs7700/LinearAlgebraUtils.jl#main",
]

foreach(url -> Pkg.add(url=url), urls)

Pkg.add([
    "ProgressMeter",
    "ProgressBars",
    "Laplacians",
    "LoggingExtras",
    "CairoMakie",
    "Makie",
])