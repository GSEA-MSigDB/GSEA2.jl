function pre_rank(fe_, sc_, se_fe_, sy_ar, ra, n_pe, n_ex, se_, ou)

    se_en = score_set(fe_, sc_, se_fe_; sy_ar...)

    _se_ra = []

    if 0 < n_pe

        println("Permuting sets to compute significance")

        se_si = Dict(se => length(fe_) for (se, fe_) in se_fe_)

        Random.seed!(ra)

        for id in ProgressBar(1:n_pe)

            se_ra = score_set(
                fe_,
                sc_,
                Dict(se => sample(fe_, si; replace = false) for (se, si) in se_si);
                sy_ar...,
            )

            push!(_se_ra, se_ra)

        end

    end

    fl_se_st = make_set_by_statistic(se_en, _se_ra, ou)

    plot_mountain(fl_se_st, n_ex, se_, fe_, sc_, se_fe_, sy_ar, ou)

    fl_se_st

end

"""
Run pre-rank GSEA

# Arguments

  - `settings_json`:
  - `set_to_genes_json`:
  - `gene_by_sample_tsv`:
  - `output_directory`:
"""
@cast function pre_rank(settings_json, set_to_genes_json, gene_by_sample_tsv, output_directory)

    ke_ar = dict_read(settings_json)

    sc_fe_sa = table_read(gene_by_sample_tsv)

    fe_ = sc_fe_sa[!, 1]

    sc_ = sc_fe_sa[!, 2]

    sc_, fe_ = sort_like([sc_, fe_])

    se_fe_ = select_set(
        dict_read(set_to_genes_json),
        ke_ar["remove_gene_set_genes"],
        fe_,
        ke_ar["minimum_gene_set_size"],
        ke_ar["maximum_gene_set_size"],
    )

    sy_ar = make_keyword_argument(ke_ar)

    pre_rank(
        fe_,
        sc_,
        se_fe_,
        sy_ar,
        ke_ar["random_seed"],
        ke_ar["number_of_permutations"],
        ke_ar["number_of_extreme_gene_sets_to_plot"],
        ke_ar["gene_sets_to_plot"],
        output_directory,
    )

end
