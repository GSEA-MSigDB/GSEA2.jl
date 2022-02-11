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

js = joinpath(TE, "set_to_genes.json")

;

GSEA.convert_gmt(joinpath(da, "h.all.v7.1.symbols.gmt"), js)

se_fe_ = OnePiece.extension.dict.read(js)

di = joinpath(da, "sarcopenia")

GSEA.convert_gct_and_cls(
    joinpath(
        di,
        "gse111016_allsamplescounts_htseqcov1_sss_forgeo.sarcopenia.vs.normal_counts_collapsed_to_symbols.gct",
    ),
    joinpath(di, "sarcopenia_binary.cls"),
    TE,
)

OnePiece.io.table.read(joinpath(TE, "score.target_by_sample.tsv"))

OnePiece.io.table.read(joinpath(TE, "score.gene_by_sample.tsv"))

GSEA.select_set(OnePiece.extension.dict.read(js), 33, 36)

se = joinpath(dirname(@__DIR__), "gsea_setting.json")

OnePiece.extension.dict.read(se)

ou = joinpath(TE, "single_sample_gsea")

mkpath(ou)

GSEA.run_single_sample_gsea(se, js, joinpath(da, "nmf_k9_w.tsv"), ou)

readdir(ou)

OnePiece.io.table.read(joinpath(ou, "enrichment.set_by_sample.tsv"))

ou = joinpath(TE, "pre_rank_gsea")

mkpath(ou)

GSEA.run_pre_rank_gsea(se, js, joinpath(da, "gene_score.tsv"), ou)

readdir(ou)

OnePiece.io.table.read(joinpath(ou, GSEA.OU))

ou = joinpath(TE, "standard_gsea")

mkpath(ou)

GSEA.run_standard_gsea(se, js, joinpath(da, "sample_value.tsv"), joinpath(da, "nmf_k9_w.tsv"), ou)

readdir(ou)

OnePiece.io.table.read(joinpath(ou, GSEA.OU))

if isdir(TE)

    rm(TE, recursive = true)

    println("Removed ", TE, ".")

end
