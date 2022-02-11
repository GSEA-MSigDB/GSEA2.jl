"""
Run single-sample GSEA

# Arguments

  - `setting_json`:
  - `set_to_genes_json`:
  - `gene_by_sample_tsv`:
  - `output_directory`:
"""
@cast function run_single_sample_gsea(setting_json, set_to_genes_json, gene_by_sample_tsv, output_directory)

    ke_ar = dict_read(setting_json)

    se_fe_ = select_set(dict_read(set_to_genes_json), pop!(ke_ar, "mi"), pop!(ke_ar, "ma"))

    sc_fe_sa = table_read(gene_by_sample_tsv)

    en_se_sa = score_set(sc_fe_sa, se_fe_; symbolize_key(ke_ar)...)

    table_write(output_directory, en_se_sa)

    en_se_sa

end
