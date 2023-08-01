using Test: @test

using BioLab

using GSEA

# ---- #

const DI = joinpath(dirname(@__DIR__), "data", "py", "1")

const IN = BioLab.Path.make_directory(joinpath(DI, "input"))

# ---- #

const TST, TSF = GSEA.convert_cls_gct(
    joinpath(IN, "target_x_sample_x_number.tsv"),
    joinpath(IN, "feature_x_sample_x_number.tsv"),
    joinpath(DI, "Coller_et_al_phen.cls"),
    joinpath(DI, "Coller_et_al_gene_exp_preproc.gct"),
)

const JS = GSEA.convert_gmt(
    joinpath(IN, "set_features.json"),
    joinpath(DI, "h.all.v2022.1.Hs.symbols.gmt"),
    joinpath(DI, "c2.all.v2022.1.Hs.symbols.gmt"),
)


# ---- #

const OU = BioLab.Path.make_directory(joinpath(DI, "output"))

# ---- #

GSEA.metric_rank(
    OU,
    TST,
    TSF,
    JS;
    normalization_dimension = 1,
    normalization_standard_deviation = 3,
    feature_x_index_x_random_tsv = joinpath(
        DI,
        "KS_SUP example mean scaling_rand_perm_gene_scores.txt",
    ),
    number_of_permutations = 1000,
)

# ---- #

const MEP = sort!(
    BioLab.DataFrame.read(joinpath(DI, "KS_SUP example mean scaling_gene_selection_scores.txt")),
    1,
)

const MEJ = sort!(BioLab.DataFrame.read(joinpath(OU, "feature_x_metric_x_score.tsv")), 1)

@test size(MEP, 1) == size(MEJ, 1)

for id in 1:size(MEP, 1)

    @test MEP[id, 1] == MEJ[id, 1]

    @test isapprox(MEP[id, 2], MEJ[id, 2]; atol = 10^-6)

end

# ---- #

const ENP = sort!(
    BioLab.DataFrame.read(joinpath(DI, "KS_SUP example mean scaling_GSEA_results_table.txt")),
    1,
)

const ENJ = sort!(BioLab.DataFrame.read(joinpath(OU, "set_x_statistic_x_number.tsv")), 1)

for id in 1:size(ENP, 1)

    @test ENP[id, 1] == ENJ[id, 1]

    @test isapprox(ENP[id, 5], ENJ[id, 2]; atol = 10^-3)

    @test isapprox(ENP[id, 4], ENJ[id, 3]; atol = 10^-2)

    @test isapprox(ENP[id, 6], ENJ[id, 4]; atol = 10^-2)

    @test isapprox(ENP[id, 7], ENJ[id, 5]; atol = 10^-2)

end
