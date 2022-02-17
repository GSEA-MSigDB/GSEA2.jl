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

println("filter!")

println("-"^99)

se_fe_ = OnePiece.dict.read(st)

co = copy(se_fe_)

GSEA.filter!(co, false, [], 33, 36)

println(keys(co))

@assert length(co) == 2

co = copy(se_fe_)

GSEA.filter!(co, true, ["SHH", "XIST"], 1, 5656)

println(co)

@assert length(co) == 2

# ----------------------------------------------------------------------------------------------- #
println("-"^99)

println("GSEA")

println("-"^99)

se = joinpath(dirname(@__DIR__), "settings.json")

println(OnePiece.dict.read(se))

# ----------------------------------------------------------------------------------------------- #
sc = joinpath(da, "score.gene_x_sample.tsv")

# ----------------------------------------------------------------------------------------------- #
println("-"^99)

println("single_sample")

println("-"^99)

ou = joinpath(TE, "single_sample_gsea")

GSEA.single_sample(se, st, sc, ou)

println(readdir(ou))

en_se_sa = OnePiece.table.read(joinpath(ou, "enrichment.set_x_sample.tsv"))

println(size(en_se_sa))

OnePiece.dataframe.view(en_se_sa)

# ----------------------------------------------------------------------------------------------- #
function print_output(ou)

    println(readdir(ou))

    fl_se_st = OnePiece.table.read(joinpath(ou, "float.set_x_statistic.tsv"))

    OnePiece.dataframe.view(fl_se_st)

end

# ----------------------------------------------------------------------------------------------- #
me = "score.gene_x_metric.tsv"

# ----------------------------------------------------------------------------------------------- #
println("-"^99)

println("pre_rank")

println("-"^99)

ou = joinpath(TE, "pre_rank_gsea")

GSEA.pre_rank(se, st, joinpath(da, me), ou)

print_output(ou)

# ----------------------------------------------------------------------------------------------- #
println("-"^99)

println("standard")

println("-"^99)

ou = joinpath(TE, "standard_gsea")

GSEA.standard(se, st, joinpath(da, "number.target_x_sample.tsv"), sc, ou)

sc_se_sa = OnePiece.table.read(joinpath(ou, me))

OnePiece.dataframe.view(sc_se_sa)

print_output(ou)

# ----------------------------------------------------------------------------------------------- #
run(`julia --project $(joinpath(@__DIR__, "sarcopenia.jl"))`)

# ----------------------------------------------------------------------------------------------- #
if isdir(TE)

    rm(TE, recursive = true)

    println("Removed $TE.")

end
