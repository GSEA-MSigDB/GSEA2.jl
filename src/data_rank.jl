"""
Run data-rank (single-sample) GSEA.

# Arguments

  - `setting_json`:
  - `gene_x_sample_x_score_tsv`:
  - `set_genes_json`:
  - `output_directory`:
"""
@cast function data_rank(setting_json, gene_x_sample_x_score_tsv, set_genes_json, output_directory)

    ke_ar = OnePiece.dict.read(setting_json)

    fe_x_sa_x_sc = OnePiece.table.read(gene_x_sample_x_score_tsv)

    se_fe_ = OnePiece.dict.read(set_genes_json)

    _filter_set!(
        se_fe_,
        ke_ar["remove_gene_set_genes"],
        fe_x_sa_x_sc[!, 1],
        ke_ar["minimum_gene_set_size"],
        ke_ar["maximum_gene_set_size"],
    )

    se_x_sa_x_en = OnePiece.feature_set_enrichment.score_set(
        fe_x_sa_x_sc,
        se_fe_;
        al = ke_ar["algorithm"],
        _make_keyword_argument(ke_ar)...,
    )

    mkpath(output_directory)

    OnePiece.table.write(joinpath(output_directory, "set_x_sample_x_enrichment.tsv"), se_x_sa_x_en)

end
