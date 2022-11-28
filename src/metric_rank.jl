function _compare_and_sort(bi_, fe_x_sa_x_sc, me, fe_)

    sc_, fes_ = BioinformaticsCore.Vector.sort_like((
        BioinformaticsCore.FeatureXSample.target(bi_, fe_x_sa_x_sc, me),
        fe_,
    ))

    fes_, sc_

end

"""
Run metric-rank (standard) GSEA.

# Arguments

  - `setting_json`:
  - `target_x_sample_x_number_tsv`:
  - `gene_x_sample_x_score_tsv`:
  - `set_genes_json`:
  - `output_directory`:
"""
@cast function metric_rank(
    setting_json,
    target_x_sample_x_number_tsv,
    gene_x_sample_x_score_tsv,
    set_genes_json,
    output_directory,
)

    #
    ke_ar = BioinformaticsCore.Dict.read(setting_json)

    #
    ta_, sat_, ta_x_sa_x_nu = BioinformaticsCore.DataFrame.separate(
        BioinformaticsCore.Table.read(target_x_sample_x_number_tsv),
    )[[2, 3, 4]]

    BioinformaticsCore.Array.error_duplicate(ta_)

    BioinformaticsCore.Matrix.error_bad(ta_x_sa_x_nu, Real)

    #
    fe_, saf_, fe_x_sa_x_sc = BioinformaticsCore.DataFrame.separate(
        BioinformaticsCore.Table.read(gene_x_sample_x_score_tsv),
    )[[2, 3, 4]]

    BioinformaticsCore.Array.error_duplicate(fe_)

    BioinformaticsCore.Matrix.error_bad(fe_x_sa_x_sc, Real)

    fe_x_sa_x_sc = fe_x_sa_x_sc[:, indexin(sat_, saf_)]

    #
    mkpath(output_directory)

    #
    bi_ = BitVector(ta_x_sa_x_nu[1, :])

    me = ke_ar["metric"]

    fe_, sc_ = _compare_and_sort(bi_, fe_x_sa_x_sc, me, fe_)

    BioinformaticsCore.Table.write(
        joinpath(output_directory, "gene_x_metric_x_score.tsv"),
        DataFrame("Gene" => fe_, me => sc_),
    )

    #
    se_fe_ = BioinformaticsCore.Dict.read(set_genes_json)

    _filter_set!(
        se_fe_,
        ke_ar["remove_gene_set_genes"],
        fe_,
        ke_ar["minimum_gene_set_size"],
        ke_ar["maximum_gene_set_size"],
    )

    #
    fe = ke_ar["feature_name"]

    sc = ke_ar["score_name"]

    al = ke_ar["algorithm"]

    sy_ar = _make_keyword_argument(ke_ar)

    pe = ke_ar["permutation"]

    ra = ke_ar["random_seed"]

    n_pe = ke_ar["number_of_permutations"]

    n_ex = ke_ar["number_of_extreme_gene_sets_to_plot"]

    pl_ = ke_ar["gene_sets_to_plot"]

    #
    if pe == "sample"

        #
        fu, id = BioinformaticsCore.FeatureSetEnrichment._match_algorithm(al)

        se_en = Dict(se => en[id] for (se, en) in fu(fe_, sc_, se_fe_; sy_ar...))

        #
        if 0 < n_pe

            println("Permuting $(pe)s to compute significance")

            #
            seed!(ra)

            #
            se_ra_ = [
                Dict(se => en[id] for (se, en) in se_en) for se_en in (
                    fu(
                        _compare_and_sort(shuffle!(bi_), fe_x_sa_x_sc, me, fe_)...,
                        se_fe_;
                        sy_ar...,
                    ) for _ in ProgressBar(1:n_pe)
                )
            ]

        else

            se_ra_ = []

        end

        #
        se_x_st_x_nu = _tabulate_statistic(se_en, se_ra_, output_directory)

        _plot_mountain(
            se_x_st_x_nu,
            fe,
            sc,
            n_ex,
            pl_,
            al,
            fe_,
            sc_,
            se_fe_,
            sy_ar,
            output_directory,
        )

        se_x_st_x_nu

        #
    elseif pe == "set"

        user_rank(fe_, sc_, se_fe_, fe, sc, al, sy_ar, ra, n_pe, n_ex, pl_, output_directory)

        #
    else

        error("`permutation` is not `sample` or `set`.")

    end

end
