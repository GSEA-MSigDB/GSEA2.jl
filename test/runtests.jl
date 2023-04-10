using Test

using BioLab

using GSEA

# --------------------------------------------- #

TE = joinpath(tempdir(), "GSEA.test")

BioLab.Path.empty(TE)

# --------------------------------------------- #

SE = joinpath(dirname(@__DIR__), "setting")

DA = @__DIR__

ST = joinpath(DA, "set_features.json")

TSS = joinpath(DA, "feature_x_sample_x_number.tsv")

TST = joinpath(DA, "target_x_sample_x_number.tsv")

TSM = "feature_x_metric_x_score.tsv"

# --------------------------------------------- #

se_fe_ = convert(Dict{String, Vector{String}}, BioLab.Dict.read(ST))

GSEA._filter_set!(se_fe_, [], 33, 36)

@test length(se_fe_) == 0

se_fe_ = convert(Dict{String, Vector{String}}, BioLab.Dict.read(ST))

GSEA._filter_set!(se_fe_, unique(vcat(values(se_fe_)...)), 33, 36)

@test length(se_fe_) == 2

se_fe_ = convert(Dict{String, Vector{String}}, BioLab.Dict.read(ST))

GSEA._filter_set!(se_fe_, ["SHH", "XIST"], 1, 5656)

@test length(se_fe_) == 2

# se_fe_ = convert(Dict{String, Vector{String}}, BioLab.Dict.read(ST))

# @code_warntype GSEA._filter_set!(se_fe_, ["SHH", "XIST"], 1, 5656)

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

    BioLab.DataFrame.print(BioLab.Table.read(joinpath(ou, "set_x_index_x_random.tsv")))

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

for (nu1_, nu2_) in (
    (collect(1:3), collect(10:10:30)),
    (fill(0, 3), fill(0, 3)),
    (fill(1, 3), fill(0.001, 3)),
    (fill(0.001, 3), fill(1, 3)),
    (fill(0.001, 3), fill(10, 3)),
)

    BioLab.print_header("$nu1_\n$nu2_")

    # TODO: `@test`.
    println(GSEA._get_signal_to_noise_ratio(nu1_, nu2_))

    # @code_warntype GSEA._get_signal_to_noise_ratio(nu1_, nu2_)

    # 22.130 ns (0 allocations: 0 bytes)
    # 10.301 ns (0 allocations: 0 bytes)
    # 28.960 ns (0 allocations: 0 bytes)
    # 28.668 ns (0 allocations: 0 bytes)
    # 28.256 ns (0 allocations: 0 bytes)
    # @btime GSEA._get_signal_to_noise_ratio($nu1_, $nu2_)

end

# --------------------------------------------- #

na = "metric_rank"

ou = joinpath(TE, na)

se = joinpath(SE, "$na.json")

GSEA.metric_rank(se, TST, TSS, ST, ou)

BioLab.DataFrame.print(BioLab.Table.read(joinpath(ou, TSM)))

print_output(ou)

readdir(joinpath(ou, "plot"))

# @code_warntype GSEA.metric_rank(se, TST, TSS, ST, ou)
