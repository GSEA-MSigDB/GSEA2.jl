using Test: @test

using BioLab

using GSEA

# ---- #

const TE = BioLab.Path.make_directory(joinpath(BioLab.TE, "GSEA"))

const DA = joinpath(dirname(@__DIR__), "data")

const PIP = joinpath(DA, "python", "input")

# ---- #

const TSF = joinpath(TE, "feature_x_sample_x_number.tsv")

const _naf, FE_, SA_, FE_X_SA_X_NU =
    BioLab.DataFrame.separate(BioLab.GCT.read(joinpath(PIP, "Coller_et_al_gene_exp_preproc.gct")))

# ---- #

foreach(BioLab.Normalization.normalize_with_0!, eachrow(FE_X_SA_X_NU))

clamp!(FE_X_SA_X_NU, -3, 3)

BioLab.DataFrame.write(TSF, BioLab.DataFrame.make("Feature", FE_, SA_, FE_X_SA_X_NU))

# ---- #

const TST = joinpath(TE, "target_x_sample_x_number.tsv")

# ---- #

const NAT, TA_, _SA_, TA_X_SA_X_NU =
    BioLab.DataFrame.separate(BioLab.CLS.read(joinpath(PIP, "Coller_et_al_phen.cls")))

BioLab.DataFrame.write(TST, BioLab.DataFrame.make(NAT, TA_, SA_, TA_X_SA_X_NU .- 1))

# ---- #

const JS = joinpath(TE, "set_features.json")

BioLab.Dict.write(
    JS,
    BioLab.Dict.merge_recursively(
        BioLab.GMT.read(joinpath(PIP, "h.all.v2022.1.Hs.symbols.gmt")),
        BioLab.GMT.read(joinpath(PIP, "c2.all.v2022.1.Hs.symbols.gmt")),
    ),
)

# ---- #

const OU = BioLab.Path.make_directory(joinpath(TE, "output"))

GSEA.metric_rank(
    OU,
    TST,
    TSF,
    JS;
    feature_x_index_x_random_tsv = joinpath(
        DA,
        "python",
        "output",
        "KS_SUP example mean scaling_rand_perm_gene_scores.txt",
    ),
    number_of_permutations = 1000,
)

# ---- #

const MEP = sort!(
    BioLab.DataFrame.read(
        joinpath(DA, "python", "output", "KS_SUP example mean scaling_gene_selection_scores.txt"),
    ),
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
    BioLab.DataFrame.read(
        joinpath(DA, "python", "output", "KS_SUP example mean scaling_GSEA_results_table.txt"),
    ),
    1,
)

ENJ = sort!(BioLab.DataFrame.read(joinpath(OU, "set_x_statistic_x_number.tsv")), 1)

for id in 1:size(ENP, 1)

    @test ENP[id, 1] == ENJ[id, 1]

    @test isapprox(ENP[id, 5], ENJ[id, 2]; atol = 10^-3)

    @test isapprox(ENP[id, 4], ENJ[id, 3]; atol = 10^-2)

    @test isapprox(ENP[id, 6], ENJ[id, 4]; atol = 10^-2)

    @test isapprox(ENP[id, 7], ENJ[id, 5]; atol = 10^-2)

end
