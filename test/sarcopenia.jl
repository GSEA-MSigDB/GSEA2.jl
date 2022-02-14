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
se = joinpath(dirname(@__DIR__), "settings.json")

ke_ar = OnePiece.extension.dict.read(se)

println(ke_ar)

# ----------------------------------------------------------------------------------------------- #
da = joinpath(@__DIR__, "sarcopenia")

# ----------------------------------------------------------------------------------------------- #
println("-"^99)

in = joinpath(da, "input")

GSEA.standard(
    se,
    joinpath(in, "set_to_genes.json"),
    joinpath(in, "int.target_by_sample.tsv"),
    joinpath(in, "float.gene_by_sample.tsv"),
    TE,
)


# ----------------------------------------------------------------------------------------------- #
be = "$(ke_ar["permutation"])_$(ke_ar["number_of_permutations"])"

# ----------------------------------------------------------------------------------------------- #
function read_old(id)

    OnePiece.io.table.read(joinpath(da, be, "gsea_report_for_$id.tsv"))

end

# ----------------------------------------------------------------------------------------------- #
ol_se_st = vcat((read_old(id) for id in [0, 1])...)

println(size(ol_se_st))

va_se_st = OnePiece.io.table.read(joinpath(TE, "score.set_by_statistic.tsv"))

println(size(va_se_st))

# ----------------------------------------------------------------------------------------------- #
se_ = ol_se_st[!, 1]

@assert isempty(symdiff(se_, va_se_st[!, 1]))

# ----------------------------------------------------------------------------------------------- #
println("-"^99)

va_se_st = va_se_st[indexin(se_, va_se_st[!, 1]), :]

for (nao, nan) in [["ES", "Enrichment"], ["NOM p-val", "P-value"], ["FDR q-val", "Q-value"]]

    display(
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
        ),
    )

end

# ----------------------------------------------------------------------------------------------- #
ol_ge_st = OnePiece.io.table.read(joinpath(da, be, "ranked_gene_list_0_versus_1.tsv"))[
    !,
    ["NAME", "SCORE"],
]

sc_ge_st = OnePiece.io.table.read(joinpath(TE, "score.gene_by_metric.tsv"))

# ----------------------------------------------------------------------------------------------- #
function ro(re)

    round(re; digits = 4)

end

# ----------------------------------------------------------------------------------------------- #
println("-"^99)

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
