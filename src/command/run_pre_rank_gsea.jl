function run_pre_rank_gsea(ke_ar, se_fe_, fe_, sc_, n_pe, ou)

    ke_ar = symbolize_key(ke_ar)

    se_en = score_set(fe_, sc_, se_fe_; ke_ar...)

    if 0 < n_pe

        println("Permuting sets to compute significance")

        se_si = Dict(se => length(fe_) for (se, fe_) in se_fe_)

        se_ra_ = []

        # TODO: set random generator

        for id in 1:n_pe

            println("  ", id, "/", n_pe)

            push!(
                se_ra_,
                score_set(
                    fe_,
                    sc_,
                    Dict(se => sample(fe_, si; replace = false) for (se, si) in se_si);
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

end

"""
Run pre-rank GSEA

# Arguments

  - `js`:
  - `se`:
  - `ts`:
  - `ou`: output directory
"""
@cast function run_pre_rank_gsea(js, se, ts, ou)

    ke_ar = dict_read(js)

    fe_sc = table_read(ts)

    run_pre_rank_gsea(
        ke_ar,
        read_set(se, ke_ar),
        fe_sc[!, 1],
        fe_sc[!, 2],
        pop!(ke_ar, "n_pe"),
        ou,
    )

end
