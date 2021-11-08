module GSEA

using Comonicon
using DataFrames
using Random
using StatsBase

using DictExtension
using FeatureBySample
using FeatureSetEnrichment
using GMTAccess
using Significance
using TableAccess

#
include("read_set.jl")

include("select_set.jl")

#
include("get_p_value_and_adjust.jl")

include("make_set_by_statistic.jl")

#
include("run_single_sample_gsea.jl")

include("run_pre_rank_gsea.jl")

include("run_standard_gsea.jl")

"""
Gene Set Enrichment Analysis.
"""
@main

end
