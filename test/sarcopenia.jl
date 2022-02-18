# ----------------------------------------------------------------------------------------------- #
TE = joinpath(tempdir(), "Sarcopenia.test")

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
da = joinpath(@__DIR__, "sarcopenia")

# ----------------------------------------------------------------------------------------------- #
println("-"^99)

println("settings.json")

println("-"^99)

se = joinpath(da, "settings.json")

ke_ar = OnePiece.dict.read(se)

println(ke_ar)

# ----------------------------------------------------------------------------------------------- #
println("-"^99)

println("GSEA standard")

println("-"^99)

GSEA.standard(
    se,
    joinpath(da, "set_genes.json"),
    joinpath(da, "number.target_x_sample.tsv"),
    joinpath(da, "score.gene_x_sample.tsv"),
    TE,
)

# ----------------------------------------------------------------------------------------------- #
be = "$(ke_ar["permutation"])_$(ke_ar["number_of_permutations"])"

# ----------------------------------------------------------------------------------------------- #
println("-"^99)

println("Set-x-Statistic")

println("-"^99)

# ----------------------------------------------------------------------------------------------- #
function re(id)

    OnePiece.table.read(joinpath(da, be, "gsea_report_for_$id.tsv"))

end

# ----------------------------------------------------------------------------------------------- #
ol_se_st = vcat((re(id) for id in [0, 1])...)

println(size(ol_se_st))

va_se_st = OnePiece.table.read(joinpath(TE, "float.set_x_statistic.tsv"))

println(size(va_se_st))

# ----------------------------------------------------------------------------------------------- #
se_ = ol_se_st[!, 1]

@assert isempty(symdiff(se_, va_se_st[!, 1]))

# ----------------------------------------------------------------------------------------------- #
println("-"^99)

va_se_st = va_se_st[indexin(se_, va_se_st[!, 1]), :]

for (nao, nan) in [
    ["ES", "Enrichment"],
    ["NES", "Normalized enrichment"],
    ["NOM p-val", "P-value"],
    #["FDR q-val", "Q-value"],
]

    OnePiece.figure.view(
        OnePiece.figure.plot_x_y(
            [ol_se_st[!, nao]],
            [va_se_st[!, nan]],
            text_ = [se_],
            mode_ = ["markers"],
            la = Dict(
                "title" => be,
                "xaxis" => Dict("title" => Dict("text" => nao)),
                "yaxis" => Dict("title" => Dict("text" => nan)),
            ),
            ou = joinpath(TE, "$(OnePiece.path.clean(nao)).html"),
        ),
    )

end

# ----------------------------------------------------------------------------------------------- #
println("-"^99)

println("Gene-x-Metric")

println("-"^99)

# ----------------------------------------------------------------------------------------------- #
ol_ge_st =
    OnePiece.table.read(joinpath(da, be, "ranked_gene_list_0_versus_1.tsv"))[!, ["NAME", "SCORE"]]

sc_ge_st = OnePiece.table.read(joinpath(TE, "score.gene_x_metric.tsv"))

# ----------------------------------------------------------------------------------------------- #
function ro(re)

    round(re; digits = 4)

end

# ----------------------------------------------------------------------------------------------- #
fe_sc = Dict(fe => ro(sc) for (fe, sc) in eachrow(sc_ge_st))

println("Gene\tOld\tNew")

n_mi = 0

for (fe, ol) in Dict(fe => ro(sc) for (fe, sc) in eachrow(ol_ge_st))

    sc = fe_sc[fe]

    if ol != sc

        println("$fe\t$ol\t$sc")

        global n_mi += 1

    end

end

n_fe = size(ol_ge_st, 1)

mi = n_mi / n_fe

println("Missed $n_mi/$n_fe scores ($(ro(mi * 100))%).")

@assert mi < 0.001

# ----------------------------------------------------------------------------------------------- #
if isdir(TE)

    rm(TE, recursive = true)

    println("Removed $TE.")

end
