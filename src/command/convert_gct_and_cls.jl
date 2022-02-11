"""
Convert `.gct` and `.cls` to `.tsv`s

# Arguments

  - `gct`:
  - `cls`:
  - `output_directory`:
"""
@cast function convert_gct_and_cls(gct, cls, output_directory)

    ge = joinpath(output_directory, "gene_by_sample.tsv")

    sc_ge_sa = gct_read(gct)

    table_write(ge, sc_ge_sa)

    ta = joinpath(output_directory, "target_by_sample.tsv")

    sa_ = names(sc_ge_sa)[2:end]

    ta_ = parse.(Float64, split(readlines(cls)[3], "\t"))

    sc_ta_sa = DataFrame(Dict(sa => ta for (sa, ta) in zip(sa_, ta_)))

    table_write(ta, sc_ta_sa)

end
