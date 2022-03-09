"""
Run single-sample GSEA

# Arguments

  - `setting_json`:
  - `set_genes_json`:
  - `gene_x_sample_x_score_tsv`:
  - `output_directory`:
"""
@cast function single_sample(
    setting_json,
    set_genes_json,
    gene_x_sample_x_score_tsv,
    output_directory,
)

    ke_ar = OnePiece.dict.read(setting_json)

    sc_fe_sa = OnePiece.table.read(gene_x_sample_x_score_tsv)

    error_feature_score(sc_fe_sa)

    se_fe_ = OnePiece.dict.read(set_genes_json)

    filter_set!(
        se_fe_,
        ke_ar["remove_gene_set_genes"],
        sc_fe_sa[!, 1],
        ke_ar["minimum_gene_set_size"],
        ke_ar["maximum_gene_set_size"],
    )

    en_se_sa = OnePiece.feature_set_enrichment.score_set(
        sc_fe_sa,
        se_fe_;
        make_keyword_argument(ke_ar)...,
    )

    mkpath(output_directory)

    OnePiece.table.write(joinpath(output_directory, "set_x_sample_x_enrichment.tsv"), en_se_sa)

end
