module GSEA

using Comonicon
using DataFrames
using OnePiece
using ProgressBars
using Random
using Statistics
using StatsBase

include("support/error_feature_score.jl")

include("support/filter!.jl")

include("support/make_keyword_argument.jl")

include("support/make_set_x_statistic.jl")

include("support/plot_mountain.jl")

include("command/single_sample.jl")

include("command/pre_rank.jl")

include("command/standard.jl")

"""
Gene Set Enrichment Analysis
"""
@main

end
