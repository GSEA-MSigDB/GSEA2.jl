function user_rank(al, fe_, sc_, se_fe_, fe, sc, sy_ar, ra, n_pe, n_ex, pl_, ou)

    se_en = BioLab.FeatureSetEnrichment.score_set(al, fe_, sc_, se_fe_; sy_ar...)

    if 0 < n_pe

        println("Permuting sets to compute significance")

        se_si = Dict(se => length(fe_) for (se, fe_) in se_fe_)

        seed!(ra)

        se_ra_ = [
            BioLab.FeatureSetEnrichment.score_set(
                al,
                fe_,
                sc_,
                Dict(se => sample(fe_, si, replace = false) for (se, si) in se_si);
                sy_ar...,
            ) for _ in ProgressBar(1:n_pe)
        ]

    else

        se_ra_ = []

    end

    se_x_st_x_nu = _tabulate_statistic(se_en, se_ra_, ou)

    _plot_mountain(se_x_st_x_nu, fe, sc, n_ex, pl_, al, fe_, sc_, se_fe_, sy_ar, ou)

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

    ke_ar = BioLab.Dict.read(setting_json)

    fe_, fe_x_me_x_sc =
        BioLab.DataFrame.separate(BioLab.Table.read(gene_x_metric_x_score_tsv))[[2, 4]]

    BioLab.Array.error_duplicate(fe_)

    BioLab.Matrix.error_bad(fe_x_me_x_sc, Real)

    sc_ = fe_x_me_x_sc[:, 1]

    sc_, fe_ = BioLab.Vector.sort_like((sc_, fe_))

    se_fe_ = BioLab.Dict.read(set_genes_json)

    _filter_set!(
        se_fe_,
        ke_ar["remove_gene_set_genes"],
        fe_,
        ke_ar["minimum_gene_set_size"],
        ke_ar["maximum_gene_set_size"],
    )

    user_rank(
        ke_ar["algorithm"],
        fe_,
        sc_,
        se_fe_,
        ke_ar["feature_name"],
        ke_ar["score_name"],
        _make_keyword_argument(ke_ar),
        ke_ar["random_seed"],
        ke_ar["number_of_permutations"],
        ke_ar["number_of_extreme_gene_sets_to_plot"],
        ke_ar["gene_sets_to_plot"],
        output_directory,
    )

end
