"""
Run standard GSEA

# Arguments

  - `setting_json`:
  - `set_to_genes_json`:
  - `target_by_sample_tsv`:
  - `gene_by_sample_tsv`:
  - `output_directory`:
"""
@cast function run_standard_gsea(
    setting_json,
    set_to_genes_json,
    target_by_sample_tsv,
    gene_by_sample_tsv,
    output_directory,
)

    ke_ar = dict_read(setting_json)

    nu_ta_sa = table_read(target_by_sample_tsv)

    sc_fe_sa = table_read(gene_by_sample_tsv)

    fe_ = string.(sc_fe_sa[:, 1])

    sc_ = compare_with_target(
        BitVector(nu_ta_sa[1, :]),
        Matrix(sc_fe_sa[:, 2:end]),
        "signal_to_noise_ratio",
    )

    se_fe_ = select_set(
        dict_read(set_to_genes_json),
        ke_ar["minimum_gene_set_size"],
        ke_ar["maximum_gene_set_size"],
    )

    pe = ke_ar["permutation"]

    n_pe = ke_ar["number_of_permutations"]

    if pe == "label"

        println("Permuting labels to compute significance")

        sy_ar = make_keyword_argument(ke_ar)

        se_en = score_set(fe_, sc_, se_fe_; sy_ar...)

        if 0 < n_pe

            sh_ = copy(sc_)

            _se_ra = []

            Random.seed!(ke_ar["random_seed"])

            for it in 1:n_pe

                println("  ", it, "/", n_pe)

                push!(_se_ra, score_set(fe_, shuffle!(sh_), se_fe_; sy_ar...))

            end

            pv_, ad_ = get_p_value_and_adjust(se_en, _se_ra)

        else

            pv_ = ad_ = fill(NaN, length(se_en))

        end

        fl_se_st = make_set_by_statistic(se_en, pv_, ad_, output_directory)

    elseif pe == "set"

        fl_se_st = run_pre_rank_gsea(ke_ar, se_fe_, fe_, sc_, n_pe, output_directory)

    else

        error("permutation is not label or set.")

    end

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
