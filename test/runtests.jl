using GSEA

using OnePiece

te = OnePiece.Path.make_temporary("GSEA.test")

se = joinpath(@__DIR__, "set_genes.json")

se_fe_ = OnePiece.Dict.read(se)

GSEA._filter_set!(se_fe_, false, [], 33, 36)

@assert length(se_fe_) == 2

se_fe_ = OnePiece.Dict.read(se)

GSEA._filter_set!(se_fe_, true, ["SHH", "XIST"], 1, 5656)

@assert length(se_fe_) == 2

GSEA._make_keyword_argument(
    Dict("exponent" => 2.0, "algorithm" => "Jensen-Shannon divergence", "number_of_jobs" => 8),
)

sc = joinpath(@__DIR__, "gene_x_sample_x_score.tsv")

ou = joinpath(te, "data_rank")

GSEA.data_rank(joinpath(dirname(@__DIR__), "setting_for_data_rank.json"), sc, se, ou)

OnePiece.DataFrame.print(OnePiece.Table.read(joinpath(ou, "set_x_sample_x_enrichment.tsv")))

me = "gene_x_metric_x_score.tsv"

function print_output(ou)

    println(readdir(ou))

    OnePiece.DataFrame.print(OnePiece.Table.read(joinpath(ou, "set_x_statistic_x_number.tsv")))

    OnePiece.DataFrame.print(OnePiece.Table.read(joinpath(ou, "set_x_random_x_enrichment.tsv")))

end

ou = joinpath(te, "user_rank")

GSEA.user_rank(
    joinpath(dirname(@__DIR__), "setting_for_user_rank.json"),
    joinpath(@__DIR__, me),
    se,
    ou,
)

print_output(ou)

readdir(joinpath(ou, "plot"))

ou = joinpath(te, "metric_rank")

GSEA.metric_rank(
    joinpath(dirname(@__DIR__), "setting_for_metric_rank.json"),
    joinpath(@__DIR__, "target_x_sample_x_number.tsv"),
    sc,
    se,
    ou,
)

OnePiece.DataFrame.print(OnePiece.Table.read(joinpath(ou, me)))

print_output(ou)

sm = joinpath(@__DIR__, "small")

GSEA.metric_rank(
    joinpath(sm, "setting_for_metric_rank.json"),
    joinpath(@__DIR__, "target_x_sample_x_number.tsv"),
    sc,
    joinpath(sm, "set_genes.json"),
    joinpath(te, "metric_rank.small"),
)
