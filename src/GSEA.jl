module GSEA

using Comonicon
using DataFrames
using OnePiece
using ProgressBars
using Random
using Statistics
using StatsBase

include("support/error_feature_score.jl")
include("support/filter_set!.jl")
include("support/make_keyword_argument.jl")
include("support/compute_statistic.jl")
include("support/plot_mountain.jl")
include("command/data_rank.jl")
include("command/user_rank.jl")
include("command/metric_rank.jl")

"""
Gene-Set Enrichment Analysis
"""
@main

end
