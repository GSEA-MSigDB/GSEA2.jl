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

TSS = joinpath(DA, "gene_x_sample_x_score.tsv")

TST = joinpath(DA, "target_x_sample_x_number.tsv")

TSM = "gene_x_metric_x_score.tsv"

# --------------------------------------------- #

se_fe_ = convert(Dict{String, Vector{String}}, BioLab.Dict.read(SG))

GSEA._filter_set!(se_fe_, false, [], 33, 36)

@test length(se_fe_) == 2

se_fe_ = convert(Dict{String, Vector{String}}, BioLab.Dict.read(SG))

GSEA._filter_set!(se_fe_, true, ["SHH", "XIST"], 1, 5656)

@test length(se_fe_) == 2

se_fe_ = convert(Dict{String, Vector{String}}, BioLab.Dict.read(SG))

# @code_warntype GSEA._filter_set!(se_fe_, true, ["SHH", "XIST"], 1, 5656)

# --------------------------------------------- #

ou = joinpath(TE, "data_rank")

GSEA.data_rank(joinpath(SE, "data_rank.json"), TSS, SG, ou)

# @code_warntype GSEA.data_rank(joinpath(SE, "data_rank.json"), TSS, SG, ou)

BioLab.DataFrame.print(BioLab.Table.read(joinpath(ou, "set_x_sample_x_enrichment.tsv")))

# --------------------------------------------- #

function print_output(ou)

    println(readdir(ou))

    BioLab.DataFrame.print(BioLab.Table.read(joinpath(ou, "set_x_statistic_x_number.tsv")))

    BioLab.DataFrame.print(BioLab.Table.read(joinpath(ou, "set_x_random_x_enrichment.tsv")))

    return nothing

end

# --------------------------------------------- #

ou = joinpath(TE, "user_rank")

GSEA.user_rank(joinpath(SE, "user_rank.json"), joinpath(DA, TSM), SG, ou)

# @code_warntype GSEA.user_rank(joinpath(SE, "user_rank.json"), joinpath(DA, TSM), SG, ou)

print_output(ou)

readdir(joinpath(ou, "plot"))

# --------------------------------------------- #

ou = joinpath(TE, "metric_rank")

GSEA.metric_rank(joinpath(SE, "metric_rank.json"), TST, TSS, SG, ou)

# @code_warntype GSEA.metric_rank(joinpath(SE, "metric_rank.json"), TST, TSS, SG, ou)

BioLab.DataFrame.print(BioLab.Table.read(joinpath(ou, TSM)))

print_output(ou)

readdir(joinpath(ou, "plot"))

sm = joinpath(DA, "small")

GSEA.metric_rank(
    joinpath(sm, "metric_rank.json"),
    TST,
    TSS,
    joinpath(sm, "set_genes.json"),
    joinpath(TE, "metric_rank.small"),
)
