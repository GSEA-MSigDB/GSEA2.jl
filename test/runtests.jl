using GSEA

using OnePiece

te = OnePiece.path.make_temporary("GSEA.test")

sett = joinpath(dirname(@__DIR__), "setting.json")

set_ = joinpath(@__DIR__, "set_genes.json")

se_fe_ = OnePiece.dict.read(set_)

GSEA._filter_set!(se_fe_, false, [], 33, 36)

@assert length(se_fe_) == 2

se_fe_ = OnePiece.dict.read(set_)

GSEA._filter_set!(se_fe_, true, ["SHH", "XIST"], 1, 5656)

@assert length(se_fe_) == 2

GSEA._make_keyword_argument(
    Dict("exponent" => 2.0, "algorithm" => "Jensen-Shannon divergence", "number_of_jobs" => 8),
)

sc = joinpath(@__DIR__, "gene_x_sample_x_score.tsv")

ou = joinpath(te, "data_rank")

GSEA.data_rank(sett, sc, set_, ou)

OnePiece.data_frame.print(OnePiece.table.read(joinpath(ou, "set_x_sample_x_enrichment.tsv")))

me = "gene_x_metric_x_score.tsv"

function print_output(ou)

    println(readdir(ou))

    OnePiece.data_frame.print(OnePiece.table.read(joinpath(ou, "set_x_statistic_x_number.tsv")))

    OnePiece.data_frame.print(OnePiece.table.read(joinpath(ou, "set_x_random_x_enrichment.tsv")))

end

OnePiece.dict.read(sett)

ou = joinpath(te, "user_rank")

GSEA.user_rank(sett, joinpath(@__DIR__, me), set_, ou)

print_output(ou)

readdir(joinpath(ou, "plot"))

ou = joinpath(te, "metric_rank")

GSEA.metric_rank(sett, joinpath(@__DIR__, "target_x_sample_x_number.tsv"), sc, set_, ou)

OnePiece.data_frame.print(OnePiece.table.read(joinpath(ou, me)))

print_output(ou)

sm = joinpath(@__DIR__, "small")

ou = joinpath(te, "metric_rank.small")

GSEA.metric_rank(
    joinpath(sm, "setting.json"),
    joinpath(@__DIR__, "target_x_sample_x_number.tsv"),
    sc,
    joinpath(sm, "set_genes.json"),
    ou,
)
