# ----------------------------------------------------------------------------------------------- #
TE = joinpath(tempdir(), "GSEA.test")

if isdir(TE)

    rm(TE, recursive = true)

    println("Removed $TE.")

end

mkdir(TE)

println("Made $TE.")

# ----------------------------------------------------------------------------------------------- #
using GSEA, OnePiece

# ----------------------------------------------------------------------------------------------- #
da = joinpath(@__DIR__, "data")

println(readdir(da))

# ----------------------------------------------------------------------------------------------- #
st = joinpath(da, "set_to_genes.json")


# ----------------------------------------------------------------------------------------------- #
println("-"^99)

se_fe_ = OnePiece.extension.dict.read(st)

in_ = ["VCAN", "SHH", "XIST"]

println(GSEA.select_set(se_fe_, false, in_, 33, 36))

println(GSEA.select_set(se_fe_, true, in_, 1, 5656))

# ----------------------------------------------------------------------------------------------- #
se = joinpath(dirname(@__DIR__), "settings.json")

println(OnePiece.extension.dict.read(se))

# ----------------------------------------------------------------------------------------------- #
fl = joinpath(da, "float.gene_by_sample.tsv")

# ----------------------------------------------------------------------------------------------- #
println("-"^99)

ou = joinpath(TE, "single_sample_gsea")

GSEA.single_sample(se, st, fl, ou)

println(readdir(ou))

println(OnePiece.io.table.read(joinpath(ou, "enrichment.set_by_sample.tsv")))

# ----------------------------------------------------------------------------------------------- #
function print_output(ou)

    println(readdir(ou))

    println(OnePiece.io.table.read(joinpath(ou, "score.set_by_statistic.tsv")))

end

# ----------------------------------------------------------------------------------------------- #
me = "score.gene_by_metric.tsv"

# ----------------------------------------------------------------------------------------------- #
println("-"^99)

ou = joinpath(TE, "pre_rank_gsea")

GSEA.pre_rank(se, st, joinpath(da, me), ou)

print_output(ou)

# ----------------------------------------------------------------------------------------------- #
println("-"^99)

ou = joinpath(TE, "standard_gsea")

GSEA.standard(se, st, joinpath(da, "int.target_by_sample.tsv"), fl, ou)

println(OnePiece.io.table.read(joinpath(ou, me)))

print_output(ou)

# ----------------------------------------------------------------------------------------------- #
if isdir(TE)

    rm(TE, recursive = true)

    println("Removed $TE.")

end
