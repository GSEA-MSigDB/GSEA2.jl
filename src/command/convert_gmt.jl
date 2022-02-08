"""
Convert `.gmt` to `.json`

# Arguments

  - `gm`: `.gmt`
  - `js`: output `.json`
"""
@cast function convert_gmt(gm, js)

    se_ge_ = gmt_read(gm)

    dict_write(js, se_ge_)

end
