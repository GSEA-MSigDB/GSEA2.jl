using Test

using BioLab

using GSEA

# --------------------------------------------- #

TE = joinpath(tempdir(), "GSEA.test")

BioLab.Path.empty(TE)

# --------------------------------------------- #

SE = joinpath(dirname(@__DIR__), "setting")

DA = @__DIR__

SG = joinpath(DA, "set_genes.json")

# --------------------------------------------- #

se_fe_ = BioLab.Dict.read(SG)

GSEA._filter_set!(se_fe_, false, [], 33, 36)

@test length(se_fe_) == 2

se_fe_ = BioLab.Dict.read(SG)

GSEA._filter_set!(se_fe_, true, ["SHH", "XIST"], 1, 5656)

@test length(se_fe_) == 2

# @code_warntype

# @btime

# --------------------------------------------- #

GSEA._make_keyword_argument(
    Dict("exponent" => 2.0, "algorithm" => "Elegant", "number_of_jobs" => 8),
)

# @code_warntype

# @btime

# --------------------------------------------- #

tss = joinpath(DA, "gene_x_sample_x_score.tsv")

# --------------------------------------------- #

ou = joinpath(TE, "data_rank")

GSEA.data_rank(joinpath(SE, "data_rank.json"), tss, SG, ou)

BioLab.DataFrame.print(BioLab.Table.read(joinpath(ou, "set_x_sample_x_enrichment.tsv")))

tsm = "gene_x_metric_x_score.tsv"

# @code_warntype

# @btime

# --------------------------------------------- #

function print_output(ou)

    println(readdir(ou))

    BioLab.DataFrame.print(BioLab.Table.read(joinpath(ou, "set_x_statistic_x_number.tsv")))

    BioLab.DataFrame.print(BioLab.Table.read(joinpath(ou, "set_x_random_x_enrichment.tsv")))

    return nothing

end

ou = joinpath(TE, "user_rank")

# --------------------------------------------- #

GSEA.user_rank(joinpath(SE, "user_rank.json"), joinpath(DA, tsm), SG, ou)

print_output(ou)

readdir(joinpath(ou, "plot"))

tst = joinpath(DA, "target_x_sample_x_number.tsv")

ou = joinpath(TE, "metric_rank")

# @code_warntype

# @btime

# --------------------------------------------- #

GSEA.metric_rank(joinpath(SE, "metric_rank.json"), tst, tss, SG, ou)

BioLab.DataFrame.print(BioLab.Table.read(joinpath(ou, tsm)))

print_output(ou)

readdir(joinpath(ou, "plot"))

sm = joinpath(DA, "small")

GSEA.metric_rank(
    joinpath(sm, "metric_rank.json"),
    tst,
    tss,
    joinpath(sm, "set_genes.json"),
    joinpath(TE, "metric_rank.small"),
)

# @code_warntype

# @btime
