module GSEA

using Comonicon
using DataFrames: DataFrame, names
using Random
using StatsBase

using OnePiece.extension.dict: read as dict_read, summarize, symbolize_key, write as dict_write
using OnePiece.extension.string: clean
using OnePiece.feature_by_sample: compare_with_target
using OnePiece.feature_set_enrichment: score_set
using OnePiece.informatics.significance: adjust_p_value, get_p_value
using OnePiece.io.gct: read as gct_read
using OnePiece.io.gmt: read as gmt_read
using OnePiece.io.table: read as table_read, write as table_write

OU = "set_by_statistic.tsv"

include("support/select_set.jl")

include("support/make_keyword_argument.jl")

include("support/get_p_value_and_adjust.jl")

include("support/make_set_by_statistic.jl")

include("support/plot_mountain.jl")

include("command/convert_gct_and_cls.jl")

include("command/convert_gmt.jl")

include("command/run_single_sample_gsea.jl")

include("command/run_pre_rank_gsea.jl")

include("command/run_standard_gsea.jl")

"""
Gene Set Enrichment Analysis
"""
@main

end
