"""
Convert `.gmt` to `.json`

# Arguments

  - `gm`: `.gmt`
  - `js`: output `.json`

"""
@cast function convert_gmt(gm::String, js::String)::Nothing

    se_ge_ = GMTAccess.read(gm)

    DictExtension.write(js, se_ge_)

    return nothing

end
