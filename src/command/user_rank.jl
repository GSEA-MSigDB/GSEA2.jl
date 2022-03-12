function user_rank(fe_, sc_, se_fe_, sy_ar, ra, n_pe, n_ex, pl_, ou)

    se_en = OnePiece.feature_set_enrichment.score_set(fe_, sc_, se_fe_; sy_ar...)

    if 0 < n_pe

        println("Permuting sets to compute significance")

        se_si = Dict(se => length(fe_) for (se, fe_) in se_fe_)

        Random.seed!(ra)

        se_ra__ = [
            OnePiece.feature_set_enrichment.score_set(
                fe_,
                sc_,
                Dict(se => sample(fe_, si, replace = false) for (se, si) in se_si);
                sy_ar...,
            ) for id in ProgressBar(1:n_pe)
        ]

    else

        se_ra__ = []

    end

    fl_se_st = compute_statistic(se_en, se_ra__, ou)

    plot_mountain(fl_se_st, n_ex, pl_, fe_, sc_, se_fe_, sy_ar, ou)

    fl_se_st

end

"""
Run user-rank (pre-rank) GSEA

# Arguments

  - `setting_json`:
  - `set_genes_json`:
  - `gene_x_metric_x_score_tsv`:
  - `output_directory`:
"""
@cast function user_rank(setting_json, set_genes_json, gene_x_metric_x_score_tsv, output_directory)

    ke_ar = OnePiece.dict.read(setting_json)

    sc_fe_sa = OnePiece.table.read(gene_x_metric_x_score_tsv)

    fe_ = sc_fe_sa[!, 1]

    sc_ = sc_fe_sa[!, 2]

    error_feature_score(fe_, sc_)

    sc_, fe_ = OnePiece.vector.sort_like(sc_, fe_)

    se_fe_ = OnePiece.dict.read(set_genes_json)

    filter_set!(
        se_fe_,
        ke_ar["remove_gene_set_genes"],
        fe_,
        ke_ar["minimum_gene_set_size"],
        ke_ar["maximum_gene_set_size"],
    )

    user_rank(
        fe_,
        sc_,
        se_fe_,
        make_keyword_argument(ke_ar),
        ke_ar["random_seed"],
        ke_ar["number_of_permutations"],
        ke_ar["number_of_extreme_gene_sets_to_plot"],
        ke_ar["gene_sets_to_plot"],
        output_directory,
    )

end
