module GSEA

using Comonicon, DataFrames, OnePiece, ProgressBars, Random, Statistics

include("support/select_set.jl")

include("support/make_keyword_argument.jl")

include("support/make_set_by_statistic.jl")

include("support/plot_mountain.jl")

include("command/single_sample.jl")

include("command/pre_rank.jl")

include("command/standard.jl")

"""
Gene Set Enrichment Analysis
"""
@main

end
