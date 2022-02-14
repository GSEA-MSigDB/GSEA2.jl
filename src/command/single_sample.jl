"""
Run single-sample GSEA

# Arguments

  - `settings_json`:
  - `set_to_genes_json`:
  - `gene_by_sample_tsv`:
  - `output_directory`:
"""
@cast function single_sample(
    settings_json,
    set_to_genes_json,
    gene_by_sample_tsv,
    output_directory,
)

    ke_ar = dict_read(settings_json)

    sc_fe_sa = table_read(gene_by_sample_tsv)

    en_se_sa = score_set(
        sc_fe_sa,
        select_set(
            dict_read(set_to_genes_json),
            ke_ar["remove_gene_set_genes"],
            sc_fe_sa[!, 1],
            ke_ar["minimum_gene_set_size"],
            ke_ar["maximum_gene_set_size"],
        );
        make_keyword_argument(ke_ar)...,
    )

    mkpath(output_directory)

    table_write(joinpath(output_directory, "enrichment.set_by_sample.tsv"), en_se_sa)

end
