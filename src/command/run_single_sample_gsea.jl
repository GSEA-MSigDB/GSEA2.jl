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

    en_se_sa = score_set(
        table_read(gene_by_sample_tsv),
        select_set(
            dict_read(set_to_genes_json),
            ke_ar["minimum_gene_set_size"],
            ke_ar["maximum_gene_set_size"],
        );
        make_keyword_argument(ke_ar)...,
    )

    table_write(joinpath(output_directory, "enrichment.set_by_sample.tsv"), en_se_sa)

end
