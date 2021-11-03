module GSEA

using Comonicon
using DictExtension
using FeatureSetEnrichment
using GMTAccess
using TableAccess

#
include("read_set.jl")

include("select.jl")

#
include("run_single_sample_gsea.jl")

include("run_pre_rank_gsea.jl")

include("run_standard_gsea.jl")

"""
Gene Set Enrichment Analysis.
"""
@main

end
