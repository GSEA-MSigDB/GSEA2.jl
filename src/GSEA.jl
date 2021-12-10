module GSEA

using Comonicon
using DataFrames
using Random
using StatsBase

using DictExtension
using FeatureBySample
using FeatureSetEnrichment
using GCTAccess
using GMTAccess
using Significance
using TableAccess

include("set/select_set.jl")

include("set/read_set.jl")

include("support/get_p_value_and_adjust.jl")

include("support/make_set_by_statistic.jl")

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
