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

    en_se_sa = score_set(table_read(gene_by_sample_tsv), se_fe_; make_keyword_argument(ke_ar)...)

    table_write(joinpath(output_directory, OU), en_se_sa)

end
