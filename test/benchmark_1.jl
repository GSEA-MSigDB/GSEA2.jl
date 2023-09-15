using Test: @test

using BioLab

using GSEA

# ---- #

const DI = joinpath(dirname(@__DIR__), "benchmark", "outer_loop")

const IN = BioLab.Path.remake_directory(joinpath(DI, "input"))

# ---- #

const TST, TSF = GSEA.convert_cls_gct(
    joinpath(IN, "target_x_sample_x_number.tsv"),
    joinpath(IN, "feature_x_sample_x_number.tsv"),
    joinpath(DI, "Coller_et_al_phen.cls"),
    joinpath(DI, "Coller_et_al_gene_exp_preproc.gct"),
)

# ---- #

const JS = GSEA.convert_gmt(
    joinpath(IN, "set_features.json"),
    joinpath(DI, "h.all.v2022.1.Hs.symbols.gmt"),
    joinpath(DI, "c2.all.v2022.1.Hs.symbols.gmt"),
)

# ---- #

const OU = GSEA.metric_rank(
    BioLab.Path.remake_directory(joinpath(DI, "output")),
    TST,
    TSF,
    JS;
    normalization_dimension = 1,
    normalization_standard_deviation = 3.0,
    permutation = joinpath(DI, "KS_SUP example mean scaling_rand_perm_gene_scores.txt"),
)

# ---- #

const PME = sort!(
    BioLab.DataFrame.read(
        joinpath(DI, "KS_SUP example mean scaling_gene_selection_scores.txt");
        select = [1, 2],
    ),
)

const JME = sort!(BioLab.DataFrame.read(joinpath(OU, "feature_x_metric_x_score.tsv")))

@test size(PME, 1) === size(JME, 1)

# ---- #

for id in 1:size(PME, 1)

    @test PME[id, 1] === JME[id, 1]

    @test isapprox(PME[id, 2], JME[id, 2]; atol = 1e-6)

end

# ---- #

const PEN = sort!(
    BioLab.DataFrame.read(
        joinpath(DI, "KS_SUP example mean scaling_GSEA_results_table.txt");
        select = [1, 5, 4, 6, 7],
    ),
)

const JEN = sort!(BioLab.DataFrame.read(joinpath(OU, "set_x_statistic_x_number.tsv")))

# ---- #

for id in 1:size(PEN, 1)

    @test PEN[id, 1] === JEN[id, 1]

    @test isapprox(PEN[id, 3], JEN[id, 2]; atol = 1e-3)

    @test isapprox(PEN[id, 2], JEN[id, 3]; atol = 1e-2)

    @test isapprox(PEN[id, 4], JEN[id, 4]; atol = 1e-2)

    @test isapprox(PEN[id, 5], JEN[id, 5]; atol = 1e-2)

end
