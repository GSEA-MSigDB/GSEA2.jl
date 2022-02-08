"""
Convert `.gct` and `.cls` to `.tsv`s

# Arguments

  - `gc`: `.gct`
  - `cl`: `.cls`
  - `ou`: output directory
"""
@cast function convert_gct_and_cls(gc::String, cl::String, ou::String)::Nothing

    sc_ge_sa = GCTAccess.read(gc)

    TableAccess.write(joinpath(ou, "score.gene_by_sample.tsv"), sc_ge_sa)

    sa_ = names(sc_ge_sa)[2:end]

    ta_ = parse.(Float64, split(readlines(cl)[3], '\t'))

    sc_ta_sa = DataFrame(Dict(sa => ta for (sa, ta) in zip(sa_, ta_)))

    TableAccess.write(joinpath(ou, "score.target_by_sample.tsv"), sc_ta_sa)::Nothing

    return nothing

end
