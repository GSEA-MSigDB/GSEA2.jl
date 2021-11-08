TE = joinpath(homedir(), "Downloads", "GSEA.test", "")

if isdir(TE)

    rm(TE; recursive = true)

end

mkdir(TE)

println("Made ", TE, ".")

using Revise
using BenchmarkTools

using GMTAccess
using TableAccess

using GSEA

da = joinpath(@__DIR__, "data", "")

gm = joinpath(da, "h.all.v7.1.symbols.gmt")

;

read_set(gm, Dict("mi" => 33, "ma" => 36))

select_set(GMTAccess.read(gm), 33, 36)

ou = joinpath(TE, "single_sample_gsea.tsv")

;

en_se_sa = run_single_sample_gsea(
    joinpath(da, "setting_for_single_sample_gsea.json"),
    gm,
    joinpath(da, "nmf_k9_w.tsv"),
    ou,
)

TableAccess.read(ou)

ou = joinpath(TE, "pre_rank_gsea", "")

mkdir(ou)

;

fl_se_st = run_pre_rank_gsea(
    joinpath(da, "setting_for_pre_rank_gsea.json"),
    gm,
    joinpath(da, "gene_score.tsv"),
    ou,
)

ou = joinpath(TE, "standard_gsea", "")

mkdir(ou)

;

fl_se_st = run_standard_gsea(
    joinpath(da, "setting_for_standard_gsea.json"),
    gm,
    joinpath(da, "sample_value.tsv"),
    joinpath(da, "nmf_k9_w.tsv"),
    ou,
)

rm(TE; recursive = true)

println("Removed ", TE, ".")
