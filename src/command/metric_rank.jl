function compare_and_sort(bi_, ma, me, fe_)

    OnePiece.vector.sort_like(OnePiece.feature_x_sample.compare_with_target(bi_, ma, me), fe_)

end

"""
Run metric-rank (standard) GSEA

# Arguments

  - `setting_json`:
  - `set_genes_json`:
  - `target_x_sample_x_number_tsv`:
  - `gene_x_sample_x_score_tsv`:
  - `output_directory`:
"""
@cast function metric_rank(
    setting_json,
    set_genes_json,
    target_x_sample_x_number_tsv,
    gene_x_sample_x_score_tsv,
    output_directory,
)

    ke_ar = OnePiece.dict.read(setting_json)

    ta_x_sa_x_nu = OnePiece.table.read(target_x_sample_x_number_tsv)

    fe_x_sa_x_sc = OnePiece.table.read(gene_x_sample_x_score_tsv)

    error_feature_score(fe_x_sa_x_sc)

    fe_ = string.(fe_x_sa_x_sc[:, 1])

    fe_x_sa_x_sc = fe_x_sa_x_sc[:, names(ta_x_sa_x_nu)]

    bi_ = BitVector(ta_x_sa_x_nu[1, :])

    ma = Matrix(fe_x_sa_x_sc)

    me = ke_ar["metric"]

    sc_, fe_ = compare_and_sort(bi_, ma, me, fe_)

    mkpath(output_directory)

    OnePiece.table.write(
        joinpath(output_directory, "gene_x_metric_x_score.tsv"),
        DataFrame("Gene" => fe_, me => sc_),
    )

    se_fe_ = OnePiece.dict.read(set_genes_json)

    filter_set!(
        se_fe_,
        ke_ar["remove_gene_set_genes"],
        fe_,
        ke_ar["minimum_gene_set_size"],
        ke_ar["maximum_gene_set_size"],
    )

    sy_ar = make_keyword_argument(ke_ar)

    pe = ke_ar["permutation"]

    ra = ke_ar["random_seed"]

    n_pe = ke_ar["number_of_permutations"]

    n_ex = ke_ar["number_of_extreme_gene_sets_to_plot"]

    pl_ = ke_ar["gene_sets_to_plot"]

    if pe == "sample"

        se_en = OnePiece.feature_set_enrichment.score_set(fe_, sc_, se_fe_; sy_ar...)

        if 0 < n_pe

            println("Permuting $(pe)s to compute significance")

            Random.seed!(ra)

            scr_, fer_ = compare_and_sort(shuffle!(bi_), ma, me, fe_)

            se_ra__ = [
                OnePiece.feature_set_enrichment.score_set(fer_, scr_, se_fe_; sy_ar...) for
                _ in ProgressBar(1:n_pe)
            ]

        else

            se_ra__ = []

        end

        se_x_st_x_nu = compute_statistic(se_en, se_ra__, output_directory)

        plot_mountain(se_x_st_x_nu, n_ex, pl_, fe_, sc_, se_fe_, sy_ar, output_directory)

        se_x_st_x_nu

    elseif pe == "set"

        user_rank(fe_, sc_, se_fe_, sy_ar, ra, n_pe, n_ex, pl_, output_directory)

    else

        error("`permutation` is not `sample` or `set`.")

    end

end
