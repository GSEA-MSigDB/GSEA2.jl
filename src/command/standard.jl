"""
Run standard GSEA

# Arguments

  - `settings_json`:
  - `set_to_genes_json`:
  - `target_by_sample_tsv`:
  - `gene_by_sample_tsv`:
  - `output_directory`:
"""
@cast function standard(
    settings_json,
    set_to_genes_json,
    target_by_sample_tsv,
    gene_by_sample_tsv,
    output_directory,
)

    ke_ar = dict_read(settings_json)

    sc_ta_sa = table_read(target_by_sample_tsv)

    sc_fe_sa = table_read(gene_by_sample_tsv)

    fe_ = string.(sc_fe_sa[:, 1])

    sc_fe_sa = sc_fe_sa[:, names(sc_ta_sa)]

    bi_ = BitVector(sc_ta_sa[1, :])

    ma = Matrix(sc_fe_sa)

    me = ke_ar["metric"]

    sc_, fe_ = sort_like([compare_with_target(bi_, ma, me), fe_])

    mkpath(output_directory)

    table_write(
        joinpath(output_directory, "gene_by_statistic.tsv"),
        DataFrame("Gene" => fe_, "Score" => sc_),
    )

    se_fe_ = select_set(
        dict_read(set_to_genes_json),
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

    se_ = ke_ar["gene_sets_to_plot"]

    if pe == "sample"

        se_en = score_set(fe_, sc_, se_fe_; sy_ar...)

        _se_ra = []

        if 0 < n_pe

            println("Permuting ", pe, "s to compute significance")

            Random.seed!(ra)

            for id in ProgressBar(1:n_pe)

                scr_, fer_ = sort_like([compare_with_target(shuffle!(bi_), ma, me), fe_])

                se_ra = score_set(fer_, scr_, se_fe_; sy_ar...)

                push!(_se_ra, se_ra)

            end

        end

        fl_se_st = make_set_by_statistic(se_en, _se_ra, output_directory)

        plot_mountain(fl_se_st, n_ex, se_, fe_, sc_, se_fe_, sy_ar, output_directory)

        fl_se_st

    elseif pe == "set"

        pre_rank(fe_, sc_, se_fe_, sy_ar, ra, n_pe, n_ex, se_, output_directory)

    else

        error("`permutation` is invalid.")

    end

end
