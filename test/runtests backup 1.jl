### A Pluto.jl notebook ###
# v0.18.0

using Markdown
using InteractiveUtils

# ╔═╡ 749f8790-8dd0-11ec-3e88-73e7bf3b61b4


# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.7.2"
manifest_format = "2.0"

[deps]
"""

# ╔═╡ Cell order:
# ╠═749f8790-8dd0-11ec-3e88-73e7bf3b61b4
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002

TE = joinpath(tempdir(), "GSEA.test")

if isdir(TE)

    rm(TE, recursive = true)

    println("Removed ", TE, ".")

end

mkdir(TE)

println("Made ", TE, ".")

using OnePiece

using GSEA

da = joinpath(@__DIR__, "data")

readdir(da)

js = joinpath(da, "set_to_genes.json")

se_ge_ = OnePiece.extension.dict.read(js)

in_ = ["VCAN", "SHH", "XIST"]

GSEA.select_set(se_ge_, false, in_, 33, 36)

GSEA.select_set(se_ge_, true, in_, 1, 5656)

se = joinpath(dirname(@__DIR__), "settings.json")

OnePiece.extension.dict.read(se)

ou = joinpath(TE, "single_sample_gsea")

GSEA.single_sample(se, js, joinpath(da, "gene_by_sample.tsv"), ou)

readdir(ou)

OnePiece.io.table.read(joinpath(ou, "enrichment.set_by_sample.tsv"))

ou = joinpath(TE, "pre_rank_gsea")

GSEA.pre_rank(se, js, joinpath(da, "gene_by_statistic.tsv"), ou)

readdir(ou)

OnePiece.io.table.read(joinpath(ou, "set_by_statistic.tsv"))

ou = joinpath(TE, "standard_gsea")

GSEA.standard(se, js, joinpath(da, "target_by_sample.tsv"), joinpath(da, "gene_by_sample.tsv"), ou)

readdir(ou)

OnePiece.io.table.read(joinpath(ou, "set_by_statistic.tsv"))

if isdir(TE)

    rm(TE, recursive = true)

    println("Removed ", TE, ".")

end
