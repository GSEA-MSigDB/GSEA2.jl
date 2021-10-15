module GSEA

include("_get_probability_and_cumulative_probability.jl")

include("make_benchmark.jl")

include("plot_mountain.jl")

include("score_set.jl")

include("score_set_new.jl")

include("try_method.jl")

using Comonicon


"""
Run single-sample GSEA.
"""
@cast function single_sample_gsea(ar; op = "Option", fl = false)

    println("ar = ", ar)

    println("op = ", op)

    println("fl = ", fl)

    return

end

"""
Run prerank GSEA.
"""
@cast function prerank_gsea(ar; op = "Option", fl = false)

    println("ar = ", ar)

    println("op = ", op)

    println("fl = ", fl)

    return

end

"""
Run GSEA.
"""
@cast function gsea(ar; op = "Option", fl = false)

    println("ar = ", ar)

    println("op = ", op)

    println("fl = ", fl)

    return

end


"""
Gene Set Enrichment Analysis.
"""
@main

end
