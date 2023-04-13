using DataFrames

using Test

using BioLab

using GSEA

# --------------------------------------------- #

di = joinpath(@__DIR__, "outer_loop")

pip = mkpath(joinpath(di, "python", "input"))

ip = mkpath(joinpath(di, "input"))

ou = mkpath(joinpath(di, "output"))

BioLab.Path.empty(ip)

BioLab.Path.empty(ou)

# --------------------------------------------- #

_fen, fe_, sa_, fe_x_sa_x_nu =
    BioLab.DataFrame.separate(BioLab.GCT.read(joinpath(pip, "Coller_et_al_gene_exp_preproc.gct")))

BioLab.Matrix.apply_by_row!(BioLab.Normalization.normalize_with_0!, fe_x_sa_x_nu)

clamp!(fe_x_sa_x_nu, -3, 3)

tsf = joinpath(ip, "feature_x_sample_x_number.tsv")

BioLab.Table.write(tsf, BioLab.DataFrame.make("Feature", fe_, sa_, fe_x_sa_x_nu))

# --------------------------------------------- #

tst = joinpath(ip, "target_x_sample_x_number.tsv")

BioLab.Table.write(
    tst,
    DataFrame(
        permutedims(
            vcat(
                "Control vs MYC",
                replace(
                    split(
                        readlines(joinpath(pip, "Coller_et_al_phen.cls"))[3],
                        ' ';
                        keepempty = false,
                    ),
                    "cntrl" => 0,
                    "myc" => 1,
                ),
            ),
        ),
        vcat("Target", sa_),
    ),
)

# --------------------------------------------- #

jss = joinpath(ip, "set_features.json")

BioLab.Dict.write(
    jss,
    BioLab.GMT.read((
        joinpath(pip, "h.all.v2022.1.Hs.symbols.gmt"),
        joinpath(pip, "c2.all.v2022.1.Hs.symbols.gmt"),
    )),
)

# --------------------------------------------- #

GSEA.metric_rank(
    joinpath(di, "setting.json"),
    tst,
    tsf,
    jss,
    ou;
    feature2_x_index_x_random = BioLab.Table.read(
        joinpath(di, "python", "txt", "KS_SUP example mean scaling_rand_perm_gene_scores.txt"),
    ),
)

# --------------------------------------------- #

dap = sort!(
    BioLab.Table.read(
        joinpath(di, "python", "txt", "KS_SUP example mean scaling_gene_selection_scores.txt"),
    ),
    1,
)

daj = sort!(BioLab.Table.read(joinpath(ou, "feature_x_metric_x_score.tsv")), 1)

@test size(dap, 1) == size(daj, 1)

for id in 1:size(dap, 1)

    @test dap[id, 1] == daj[id, 1]

    @test isapprox(dap[id, 2], daj[id, 2]; atol = 10^-6)

end

# --------------------------------------------- #

dap = sort!(
    BioLab.Table.read(
        joinpath(di, "python", "txt", "KS_SUP example mean scaling_GSEA_results_table.txt"),
    ),
    1,
)

daj = sort!(BioLab.Table.read(joinpath(ou, "set_x_statistic_x_number.tsv")), 1)

for id in 1:size(dap, 1)

    @test dap[id, 1] == daj[id, 1]

    @test isapprox(dap[id, 5], daj[id, 2]; atol = 10^-3)

    @test isapprox(dap[id, 4], daj[id, 3]; atol = 10^-2)

    @test isapprox(dap[id, 6], daj[id, 4]; atol = 10^-2)

    @test isapprox(dap[id, 7], daj[id, 5]; atol = 10^-2)

end
