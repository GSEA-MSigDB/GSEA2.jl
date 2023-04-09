using Test

using BioLab

using GSEA

# --------------------------------------------- #

TE = joinpath(tempdir(), "GSEA.test")

BioLab.Path.empty(TE)

# --------------------------------------------- #

SE = joinpath(dirname(@__DIR__), "setting")

DA = @__DIR__

ST = joinpath(DA, "set_genes.json")

TSS = joinpath(DA, "gene_x_sample_x_score.tsv")

TST = joinpath(DA, "target_x_sample_x_number.tsv")

TSM = "gene_x_metric_x_score.tsv"

# --------------------------------------------- #

se_fe_ = convert(Dict{String, Vector{String}}, BioLab.Dict.read(ST))

GSEA._filter_set!(se_fe_, false, [], 33, 36)

@test length(se_fe_) == 2

se_fe_ = convert(Dict{String, Vector{String}}, BioLab.Dict.read(ST))

GSEA._filter_set!(se_fe_, true, ["SHH", "XIST"], 1, 5656)

@test length(se_fe_) == 2

# se_fe_ = convert(Dict{String, Vector{String}}, BioLab.Dict.read(ST))

# @code_warntype GSEA._filter_set!(se_fe_, true, ["SHH", "XIST"], 1, 5656)

# --------------------------------------------- #

na = "data_rank"

ou = joinpath(TE, na)

se = joinpath(SE, "$na.json")

GSEA.data_rank(se, TSS, ST, ou)

BioLab.DataFrame.print(BioLab.Table.read(joinpath(ou, "set_x_sample_x_enrichment.tsv")))

# @code_warntype GSEA.data_rank(se, TSS, ST, ou)

# --------------------------------------------- #

function print_output(ou)

    println(readdir(ou))

    BioLab.DataFrame.print(BioLab.Table.read(joinpath(ou, "set_x_statistic_x_number.tsv")))

    BioLab.DataFrame.print(BioLab.Table.read(joinpath(ou, "set_x_random_x_enrichment.tsv")))

    return nothing

end

# --------------------------------------------- #

na = "user_rank"

ou = joinpath(TE, na)

se = joinpath(SE, "$na.json")

me = joinpath(DA, TSM)

GSEA.user_rank(se, me, ST, ou)

print_output(ou)

readdir(joinpath(ou, "plot"))

# @code_warntype GSEA.user_rank(se, me, ST, ou)

# --------------------------------------------- #

na = "metric_rank"

ou = joinpath(TE, na)

se = joinpath(SE, "$na.json")

GSEA.metric_rank(se, TST, TSS, ST, ou)

BioLab.DataFrame.print(BioLab.Table.read(joinpath(ou, TSM)))

print_output(ou)

readdir(joinpath(ou, "plot"))

# @code_warntype GSEA.metric_rank(se, TST, TSS, ST, ou)

# --------------------------------------------- #

sm = joinpath(DA, "small")

GSEA.metric_rank(
    joinpath(sm, "setting.json"),
    TST,
    TSS,
    joinpath(sm, "set_genes.json"),
    joinpath(TE, "metric_rank.small"),
)
