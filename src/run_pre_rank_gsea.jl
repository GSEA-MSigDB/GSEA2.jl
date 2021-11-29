function run_pre_rank_gsea(
    ke_ar::Dict{String, Any},
    se_fe_::Dict{String, Vector{String}},
    fe_::Vector{String},
    sc_::Vector{Float64},
    n_pe::Int64,
    ou::String,
)::DataFrame

    ke_ar = DictExtension.convert_to_keyword_argument(ke_ar)

    se_en = FeatureSetEnrichment.score_set(fe_, sc_, se_fe_; ke_ar...)

    if 0 < n_pe

        println("Permuting sets to compute significance")

        se_si = Dict(se => length(fe_) for (se, fe_) in se_fe_)

        se_ra_ = Vector{Dict{String, Float64}}()

        println("SET UP RANDOM GENERATOR")

        for id in 1:n_pe

            println("  ", id, "/", n_pe)

            push!(
                se_ra_,
                FeatureSetEnrichment.score_set(
                    fe_,
                    sc_,
                    Dict(
                        se => sample(fe_, si; replace = false) for
                        (se, si) in se_si
                    );
                    ke_ar...,
                ),
            )

        end

        pv_, ad_ = get_p_value_and_adjust(se_en, se_ra_)

    else

        pv_ = ad_ = fill(NaN, length(se_en))

    end

    fl_se_st = make_set_by_statistic(se_en, pv_, ad_, ou)

    println("Plotting")

    return fl_se_st


end

"""
Run pre-rank GSEA
"""
@cast function run_pre_rank_gsea(
    js::String,
    gm::String,
    ts::String,
    ou::String,
)::DataFrame

    ke_ar = DictExtension.read(js)

    fe_sc = TableAccess.read(ts)

    return run_pre_rank_gsea(
        ke_ar,
        read_set(gm, ke_ar),
        convert(Vector{String}, fe_sc[!, 1]),
        fe_sc[!, 2],
        pop!(ke_ar, "n_pe"),
        ou,
    )

end
