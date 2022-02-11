"""
Convert `.gmt` to `.json`

# Arguments

  - `gmt`:
  - `set_to_genes_json`:
"""
@cast function convert_gmt(gmt, set_to_genes_json)

    se_ge_ = gmt_read(gmt)

    dict_write(set_to_genes_json, se_ge_)

end
