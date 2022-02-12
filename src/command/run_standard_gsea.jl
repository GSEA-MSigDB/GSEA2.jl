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

    sc_ta_sa = table_read(target_by_sample_tsv)

    sc_fe_sa = table_read(gene_by_sample_tsv)

    fe_ = string.(sc_fe_sa[:, 1])

    sc_ = compare_with_target(
        BitVector(sc_ta_sa[1, :]),
        Matrix(sc_fe_sa[:, 2:end]),
        "signal_to_noise_ratio",
    )

    sc_, fe_ = sort_like([sc_, fe_])

    se_fe_ = select_set(
        dict_read(set_to_genes_json),
        ke_ar["minimum_gene_set_size"],
        ke_ar["maximum_gene_set_size"],
    )

    sy_ar = make_keyword_argument(ke_ar)

    pe = ke_ar["permutation"]

    ra = ke_ar["random_seed"]

    n_pe = ke_ar["number_of_permutations"]

    if pe == "label"

        se_en = score_set(fe_, sc_, se_fe_; sy_ar...)

        ra__ = []

        if 0 < n_pe

            println("Permuting labels to compute significance")

            sh_ = copy(sc_)

            Random.seed!(ra)

            pr = round(n_pe / 10)

            for id in 1:n_pe

                if convert(Bool, id % pr)

                    println("  ", id, "/", n_pe)

                end

                push!(ra__, collect(values(score_set(fe_, shuffle!(sh_), se_fe_; sy_ar...))))

            end

        end

        fl_se_st = make_set_by_statistic(se_en, ra__, output_directory)

        plot_mountain(
            fl_se_st,
            ke_ar["number_of_extreme_gene_sets_to_plot"],
            ke_ar["gene_sets_to_plot"],
            fe_,
            sc_,
            se_fe_,
            sy_ar,
            output_directory,
        )

        fl_se_st

    elseif pe == "set"

        run_pre_rank_gsea(fe_, sc_, se_fe_, sy_ar, ra, n_pe, output_directory)

    else

        error("permutation is not label or set.")

    end

end
