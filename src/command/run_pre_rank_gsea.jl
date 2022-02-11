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

    for ro in eachrow(fl_se_st)

        println(ro[1])

    end

end

"""
Run pre-rank GSEA

# Arguments

  - `setting_json`:
  - `set_to_genes_json`:
  - `gene_by_sample_tsv`:
  - `output_directory`:
"""
@cast function run_pre_rank_gsea(
    setting_json,
    set_to_genes_json,
    gene_by_sample_tsv,
    output_directory,
)

    ke_ar = dict_read(setting_json)

    se_fe_ = dict_read(set_to_genes_json)

    se_fe_ = select_set(se_fe_, pop!(ke_ar, "mi"), pop!(ke_ar, "ma"))

    fe_sc = table_read(gene_by_sample_tsv)

    run_pre_rank_gsea(
        ke_ar,
        se_fe_,
        fe_sc[!, 1],
        fe_sc[!, 2],
        pop!(ke_ar, "n_pe"),
        output_directory,
    )

end
