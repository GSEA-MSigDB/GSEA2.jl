module GSEA

using Comonicon
using DataFrames: DataFrame, names
using ProgressBars
using Random
using StatsBase

using OnePiece.extension.dict: read as dict_read, summarize, symbolize_key, write as dict_write
using OnePiece.extension.path: clean
using OnePiece.extension.vector: sort_like
using OnePiece.feature_by_sample: compare_with_target
using OnePiece.feature_set_enrichment: score_set
using OnePiece.informatics.significance: adjust_p_value, get_p_value
using OnePiece.io.table: read as table_read, write as table_write

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
