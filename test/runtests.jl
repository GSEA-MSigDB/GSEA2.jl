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

se = joinpath(dirname(@__DIR__), "settings.json")

OnePiece.extension.dict.read(se)

ou = joinpath(TE, "single_sample_gsea")

mkpath(ou)

GSEA.single_sample(se, js, joinpath(da, "gene_by_sample.tsv"), ou)

readdir(ou)

OnePiece.io.table.read(joinpath(ou, "enrichment.set_by_sample.tsv"))

ou = joinpath(TE, "pre_rank_gsea")

mkpath(ou)

GSEA.pre_rank(se, js, joinpath(da, "gene_by_statistic.tsv"), ou)

readdir(ou)

OnePiece.io.table.read(joinpath(ou, "set_by_statistic.tsv"))

ou = joinpath(TE, "standard_gsea")

mkpath(ou)

GSEA.standard(se, js, joinpath(da, "target_by_sample.tsv"), joinpath(da, "gene_by_sample.tsv"), ou)

readdir(ou)

OnePiece.io.table.read(joinpath(ou, "set_by_statistic.tsv"))

da = joinpath(@__DIR__, "sarcopenia")

target_by_sample_tsv = joinpath(da, "target_by_sample.tsv")

gene_by_sample_tsv = joinpath(da, "gene_by_sample.tsv")

readdir(da)

ou = joinpath(TE, "sarcopenia")

mkpath(ou)

GSEA.standard(se, joinpath(da, "set_to_genes.json"), target_by_sample_tsv, gene_by_sample_tsv, ou)

readdir(ou)

na = "gene_by_statistic.tsv"

ol_ge_st = OnePiece.io.table.read(joinpath(da, na))

nu_ge_st = OnePiece.io.table.read(joinpath(ou, na))

fe_sc = Dict(fe => sc for (fe, sc) in eachrow(nu_ge_st))

println("Gene", "\t", "Old", "\t", "New")

di = 4

n_mi = 0

for (id, (fe, ol)) in enumerate(eachrow(ol_ge_st))

    ne = fe_sc[fe]

    ol = round(ol; digits = di)

    ne = round(ne; digits = di)

    if ol != ne

        println(fe, "\t", ol, "\t", ne)

        n_mi += 1

    end

end

n_fe = size(ol_ge_st)[1]

mi = n_mi / n_fe

println("Missed ", n_mi, "/", n_fe, " scores (", round(mi * 100; digits = di), "%).")

@assert mi < 0.001

nu_se_st = OnePiece.io.table.read(joinpath(ou, "set_by_statistic.tsv"))

if isdir(TE)

    rm(TE, recursive = true)

    println("Removed ", TE, ".")

end
