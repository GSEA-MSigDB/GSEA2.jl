TE = joinpath(tempdir(), "GSEA.test")

if isdir(TE)

    rm(TE; recursive = true)

end

mkdir(TE)

#using Revise
#using BenchmarkTools

using OnePiece

using GSEA

da = joinpath(@__DIR__, "data")

se = joinpath(da, "h.all.v7.1.symbols.gmt")

;

GSEA.select_set(OnePiece.io.gmt.read(se), 33, 36)

GSEA.read_set(se, Dict("mi" => 33, "ma" => 36))

ou = joinpath(TE, "single_sample_gsea.tsv")

GSEA.run_single_sample_gsea(
    joinpath(da, "setting", "single_sample_gsea.json"),
    se,
    joinpath(da, "nmf_k9_w.tsv"),
    ou,
)

readdir(TE)

ou = joinpath(TE, "pre_rank_gsea")

mkdir(ou)

GSEA.run_pre_rank_gsea(
    joinpath(da, "setting", "pre_rank_gsea.json"),
    se,
    joinpath(da, "gene_score.tsv"),
    ou,
)

readdir(ou)

ou = joinpath(TE, "standard_gsea")

mkdir(ou)

GSEA.run_standard_gsea(
    joinpath(da, "setting", "standard_gsea.json"),
    se,
    joinpath(da, "sample_value.tsv"),
    joinpath(da, "nmf_k9_w.tsv"),
    ou,
)

readdir(ou)

di = joinpath(da, "sarcopenia")

gc = joinpath(
    di,
    "gse111016_allsamplescounts_htseqcov1_sss_forgeo.sarcopenia.vs.normal_counts_collapsed_to_symbols.gct",
)

cl = joinpath(di, "sarcopenia_binary.cls")

GSEA.convert_gct_and_cls(gc, cl, di)

readdir(di)

gm = joinpath(di, "c2.cp.wikipathways.v7.4.symbols.gmt")

js = joinpath(di, "set_to_genes.json")

GSEA.convert_gmt(gm, js)

readdir(di)

rm(TE; recursive = true)
