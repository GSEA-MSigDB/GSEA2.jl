using GSEA

using BioLab

se = joinpath(dirname(@__DIR__), "setting")

da = joinpath(@__DIR__, "data")

js = joinpath(da, "set_genes.json")

te = BioLab.Path.make_temporary("GSEA.test")

;

se_fe_ = BioLab.Dict.read(js)

GSEA._filter_set!(se_fe_, false, [], 33, 36)

if length(se_fe_) != 2

    error()

end

se_fe_ = BioLab.Dict.read(js)

GSEA._filter_set!(se_fe_, true, ["SHH", "XIST"], 1, 5656)

if length(se_fe_) != 2

    error()

end

GSEA._make_keyword_argument(
    Dict("exponent" => 2.0, "algorithm" => "Jensen-Shannon divergence", "number_of_jobs" => 8),
)

tss = joinpath(da, "gene_x_sample_x_score.tsv")

ou = joinpath(te, "data_rank")

;

GSEA.data_rank(joinpath(se, "data_rank.json"), tss, js, ou)

BioLab.DataFrame.print(BioLab.Table.read(joinpath(ou, "set_x_sample_x_enrichment.tsv")))

tsm = "gene_x_metric_x_score.tsv"

;

function print_output(ou)

    println(readdir(ou))

    BioLab.DataFrame.print(BioLab.Table.read(joinpath(ou, "set_x_statistic_x_number.tsv")))

    BioLab.DataFrame.print(BioLab.Table.read(joinpath(ou, "set_x_random_x_enrichment.tsv")))

end

ou = joinpath(te, "user_rank")

;

GSEA.user_rank(joinpath(se, "user_rank.json"), joinpath(da, tsm), js, ou)

print_output(ou)

readdir(joinpath(ou, "plot"))

tst = joinpath(da, "target_x_sample_x_number.tsv")

ou = joinpath(te, "metric_rank")

;

GSEA.metric_rank(joinpath(se, "metric_rank.json"), tst, tss, js, ou)

BioLab.DataFrame.print(BioLab.Table.read(joinpath(ou, tsm)))

print_output(ou)

readdir(joinpath(ou, "plot"))

sm = joinpath(da, "small")

GSEA.metric_rank(
    joinpath(sm, "metric_rank.json"),
    tst,
    tss,
    joinpath(sm, "set_genes.json"),
    joinpath(te, "metric_rank.small"),
)
