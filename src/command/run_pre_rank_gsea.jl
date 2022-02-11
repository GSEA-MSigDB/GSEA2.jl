function run_pre_rank_gsea(ke_ar, se_fe_, fe_, sc_, n_pe, ou)

    sy_ar = make_keyword_argument(ke_ar)

    se_en = score_set(fe_, sc_, se_fe_; sy_ar...)

    if 0 < n_pe

        println("Permuting sets to compute significance")

        se_si = Dict(se => length(fe_) for (se, fe_) in se_fe_)

        se_ra_ = []

        Random.seed!(ke_ar["random_seed"])

        for id in 1:n_pe

            println("  ", id, "/", n_pe)

            push!(
                se_ra_,
                score_set(
                    fe_,
                    sc_,
                    Dict(se => sample(fe_, si; replace = false) for (se, si) in se_si);
                    sy_ar...,
                ),
            )

        end

        pv_, ad_ = get_p_value_and_adjust(se_en, se_ra_)

    else

        pv_ = ad_ = fill(NaN, length(se_en))

    end

    make_set_by_statistic(se_en, pv_, ad_, ou)

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

    fe_sc = table_read(gene_by_sample_tsv)

    fl_se_st = run_pre_rank_gsea(
        ke_ar,
        select_set(
            dict_read(set_to_genes_json),
            ke_ar["minimum_gene_set_size"],
            ke_ar["maximum_gene_set_size"],
        ),
        fe_sc[!, 1],
        fe_sc[!, 2],
        ke_ar["number_of_permutations"],
        output_directory,
    )

    plot_mountain(
        fl_se_st,
        ke_ar["number_of_extreme_gene_sets_to_plot"],
        ke_ar["gene_sets_to_plot"],
        sc_fe_sa,
        se_fe_,
        sy_ar,
        output_directory,
    )

end
