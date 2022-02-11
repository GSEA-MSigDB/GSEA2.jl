"""
Run single-sample GSEA

# Arguments

  - `setting_json`:
  - `set_to_genes_json`:
  - `gene_by_sample_tsv`:
  - `output_directory`:
"""
@cast function run_single_sample_gsea(
    setting_json,
    set_to_genes_json,
    gene_by_sample_tsv,
    output_directory,
)

    ke_ar = dict_read(setting_json)

    se_fe_ = select_set(
        dict_read(set_to_genes_json),
        ke_ar["minimum_gene_set_size"],
        ke_ar["maximum_gene_set_size"],
    )

    sy_ar = make_keyword_argument(ke_ar)

    sc_fe_sa = table_read(gene_by_sample_tsv)

    en_se_sa = score_set(sc_fe_sa, se_fe_; sy_ar...)

    table_write(joinpath(output_directory, OU), en_se_sa)

    plot_mountain(
        en_se_sa,
        ke_ar["number_of_extreme_gene_sets_to_plot"],
        ke_ar["gene_sets_to_plot"],
        sc_fe_sa,
        se_fe_,
        sy_ar,
    )

end
