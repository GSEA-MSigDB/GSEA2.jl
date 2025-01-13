module GSEA

using ProgressMeter: @showprogress

using Omics

include("enrich.jl")

include("plot.jl")

include("rank.jl")

end
