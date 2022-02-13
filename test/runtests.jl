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

GSEA.select_set(OnePiece.extension.dict.read(js), 33, 36)

GSEA.select_set(OnePiece.extension.dict.read(js), 1, 5656, in_ = ["VCAN", "SHH", "XIST"])

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
