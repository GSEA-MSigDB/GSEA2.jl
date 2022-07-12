module GSEA

using Comonicon
using DataFrames
using OnePiece
using ProgressBars
using Random
using Statistics
using StatsBase

include("_compute_statistic.jl")

include("_error_feature_score.jl")

include("_filter_set!.jl")

include("_make_keyword_argument.jl")

include("_plot_mountain.jl")

include("data_rank.jl")

include("metric_rank.jl")

include("user_rank.jl")

"""
Gene-Set Enrichment Analysis
"""
@main

end
