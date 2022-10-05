function user_rank(fe_, sc_, se_fe_, al, sy_ar, ra, n_pe, n_ex, pl_, ou)

    fu, st = OnePiece.feature_set_enrichment._match_algorithm(al)

    se_en = OnePiece.feature_set_enrichment._match_algorithm(fu(fe_, sc_, se_fe_; sy_ar...), st)

    if 0 < n_pe

        println("Permuting sets to compute significance")

        se_le = Dict(se => length(fe_) for (se, fe_) in se_fe_)

        Random.seed!(ra)

        se_ra__ = OnePiece.feature_set_enrichment._match_algorithm(
            [
                fu(
                    fe_,
                    sc_,
                    Dict(se => sample(fe_, si, replace = false) for (se, si) in se_le);
                    sy_ar...,
                ) for _ in ProgressBar(1:n_pe)
            ],
            st,
        )

    else

        se_ra__ = []

    end

    se_x_st_x_nu = _compute_statistic(se_en, se_ra__, ou)

    _plot_mountain(se_x_st_x_nu, n_ex, pl_, al, fe_, sc_, se_fe_, sy_ar, ou)

    se_x_st_x_nu

end

"""
Run user-rank (pre-rank) GSEA.

# Arguments

  - `setting_json`:
  - `gene_x_metric_x_score_tsv`:
  - `set_genes_json`:
  - `output_directory`:
"""
@cast function user_rank(setting_json, gene_x_metric_x_score_tsv, set_genes_json, output_directory)

    ke_ar = OnePiece.dict.read(setting_json)

    fe_, ma = OnePiece.data_frame.separate(OnePiece.table.read(gene_x_metric_x_score_tsv))[[2, 4]]

    OnePiece.vector.error_duplicate(fe_)

    OnePiece.matrix.error_bad(ma, Real)

    sc_ = ma[:, 1]

    sc_, fe_ = OnePiece.vector.sort_like(sc_, fe_)

    se_fe_ = OnePiece.dict.read(set_genes_json)

    _filter_set!(
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
        ke_ar["algorithm"],
        _make_keyword_argument(ke_ar),
        ke_ar["random_seed"],
        ke_ar["number_of_permutations"],
        ke_ar["number_of_extreme_gene_sets_to_plot"],
        ke_ar["gene_sets_to_plot"],
        output_directory,
    )

end
