# ----------------------------------------------------------------------------------------------- #
TE = joinpath(tempdir(), "GSEA.test")

if isdir(TE)

    rm(TE, recursive = true)

    println("Removed $TE.")

end

mkdir(TE)

println("Made $TE.")

# ----------------------------------------------------------------------------------------------- #
using GSEA
using OnePiece

# ----------------------------------------------------------------------------------------------- #
da = joinpath(@__DIR__, "data")

st = joinpath(da, "set_genes.json")

# ----------------------------------------------------------------------------------------------- #
println("-"^99)

println("filter_set!")

println("-"^99)

se_fe_ = OnePiece.dict.read(st)

co = copy(se_fe_)

GSEA.filter_set!(co, false, [], 33, 36)

println(keys(co))

@assert length(co) == 2

co = copy(se_fe_)

GSEA.filter_set!(co, true, ["SHH", "XIST"], 1, 5656)

println(co)

@assert length(co) == 2

# ----------------------------------------------------------------------------------------------- #
println("-"^99)

println("make_keyword_argument")

println("-"^99)

println(
    GSEA.make_keyword_argument(
        Dict("exponent" => 2.0, "algorithm" => "Jensen-Shannon divergence", "number_of_jobs" => 8),
    ),
)

# ----------------------------------------------------------------------------------------------- #
println("-"^99)

println("GSEA")

println("-"^99)

se = joinpath(dirname(@__DIR__), "setting.json")

println(OnePiece.dict.read(se))

# ----------------------------------------------------------------------------------------------- #
sc = joinpath(da, "gene_x_sample_x_score.tsv")

# ----------------------------------------------------------------------------------------------- #
su = "data_rank"

println("-"^99)

println(su)

println("-"^99)

ou = joinpath(TE, su)

GSEA.data_rank(se, st, sc, ou)

println(readdir(ou))

en_se_sa = OnePiece.table.read(joinpath(ou, "set_x_sample_x_enrichment.tsv"))

println(size(en_se_sa))

OnePiece.dataframe.view(en_se_sa)

# ----------------------------------------------------------------------------------------------- #
function print_output(ou)

    println(readdir(ou))

    fl_se_st = OnePiece.table.read(joinpath(ou, "set_x_statistic_x_number.tsv"))

    OnePiece.dataframe.view(fl_se_st)

end

# ----------------------------------------------------------------------------------------------- #
me = "gene_x_metric_x_score.tsv"

# ----------------------------------------------------------------------------------------------- #
su = "user_rank"

println("-"^99)

println(su)

println("-"^99)

ou = joinpath(TE, su)

GSEA.user_rank(se, st, joinpath(da, me), ou)

print_output(ou)

# ----------------------------------------------------------------------------------------------- #
su = "metric_rank"

println("-"^99)

println(su)

println("-"^99)

ou = joinpath(TE, su)

GSEA.metric_rank(se, st, joinpath(da, "target_x_sample_x_number.tsv"), sc, ou)

sc_se_sa = OnePiece.table.read(joinpath(ou, me))

OnePiece.dataframe.view(sc_se_sa)

print_output(ou)

# ----------------------------------------------------------------------------------------------- #
run(`julia --project $(joinpath(@__DIR__, "sarcopenia.jl"))`)

# ----------------------------------------------------------------------------------------------- #
println("Done.")
