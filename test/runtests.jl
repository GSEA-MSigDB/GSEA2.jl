using GSEA
using OnePiece

TE = OnePiece.path.make_temporary("GSEA.test")

sett = joinpath(dirname(@__DIR__), "setting.json")

set_ = joinpath(@__DIR__, "set_genes.json")

te = joinpath(tempdir(), "GSEA")

rm(te, recursive = true)

mkdir(te)

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

GSEA.data_rank(sett, set_, sc, ou)

OnePiece.dataframe.print(OnePiece.table.read(joinpath(ou, "set_x_sample_x_enrichment.tsv")))

me = "gene_x_metric_x_score.tsv"

function print_output(ou)

    println(readdir(ou))

    OnePiece.dataframe.print(OnePiece.table.read(joinpath(ou, "set_x_statistic_x_number.tsv")))

    OnePiece.dataframe.print(OnePiece.table.read(joinpath(ou, "set_x_random_x_enrichment.tsv")))

end

ou = joinpath(te, "user_rank")

GSEA.user_rank(sett, set_, joinpath(@__DIR__, me), ou)

print_output(ou)

ou = joinpath(te, "metric_rank")

GSEA.metric_rank(sett, set_, joinpath(@__DIR__, "target_x_sample_x_number.tsv"), sc, ou)

OnePiece.dataframe.print(OnePiece.table.read(joinpath(ou, me)))

print_output(ou)
