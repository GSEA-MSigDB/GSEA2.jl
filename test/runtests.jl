using Test

using BioLab

using GSEA

# ----------------------------------------------------------------------------------------------- #

se = joinpath(dirname(@__DIR__), "setting")

da = joinpath(@__DIR__, "data")

sg = joinpath(da, "set_genes.json")

te = BioLab.Path.make_temporary("GSEA.test")

# ----------------------------------------------------------------------------------------------- #

se_fe_ = BioLab.Dict.read(sg)

GSEA._filter_set!(se_fe_, false, [], 33, 36)

@test length(se_fe_) == 2

se_fe_ = BioLab.Dict.read(sg)

GSEA._filter_set!(se_fe_, true, ["SHH", "XIST"], 1, 5656)

@test length(se_fe_) == 2

# @code_warntype

# @btime

# ----------------------------------------------------------------------------------------------- #

GSEA._make_keyword_argument(
    Dict("exponent" => 2.0, "algorithm" => "Elegant", "number_of_jobs" => 8),
)

# @code_warntype

# @btime

# ----------------------------------------------------------------------------------------------- #

tss = joinpath(da, "gene_x_sample_x_score.tsv")

# ----------------------------------------------------------------------------------------------- #

ou = joinpath(te, "data_rank")

GSEA.data_rank(joinpath(se, "data_rank.json"), tss, sg, ou)

BioLab.DataFrame.print(BioLab.Table.read(joinpath(ou, "set_x_sample_x_enrichment.tsv")))

tsm = "gene_x_metric_x_score.tsv"

# @code_warntype

# @btime

# ----------------------------------------------------------------------------------------------- #

function print_output(ou)

    println(readdir(ou))

    BioLab.DataFrame.print(BioLab.Table.read(joinpath(ou, "set_x_statistic_x_number.tsv")))

    BioLab.DataFrame.print(BioLab.Table.read(joinpath(ou, "set_x_random_x_enrichment.tsv")))

    return nothing

end

ou = joinpath(te, "user_rank")

# ----------------------------------------------------------------------------------------------- #

GSEA.user_rank(joinpath(se, "user_rank.json"), joinpath(da, tsm), sg, ou)

print_output(ou)

readdir(joinpath(ou, "plot"))

tst = joinpath(da, "target_x_sample_x_number.tsv")

ou = joinpath(te, "metric_rank")

# @code_warntype

# @btime

# ----------------------------------------------------------------------------------------------- #

GSEA.metric_rank(joinpath(se, "metric_rank.json"), tst, tss, sg, ou)

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

# @code_warntype

# @btime
