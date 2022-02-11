"""
Convert `.gmt` to `.json`

# Arguments

  - `gmt`:
  - `set_to_genes_json`:
"""
@cast function convert_gmt(gmt, set_to_genes_json)

    dict_write(set_to_genes_json, gmt_read(gmt))

end
