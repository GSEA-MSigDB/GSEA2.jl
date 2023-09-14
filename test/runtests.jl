using Aqua: test_all, test_ambiguities

using Test: @test

using GSEA

# ---- #

test_all(GSEA; ambiguities = false)

test_ambiguities(GSEA)

# ----------------------------------------------------------------------------------------------- #

#run(`julia --project FeatureSetEnrichment.jl`)

# ---- #

using BioLab

# ---- #

const TE = BioLab.Path.remake_directory(joinpath(BioLab.TE, "GSEA"))

# ---- #

const DA = joinpath(dirname(@__DIR__), "data")

@test BioLab.Path.read(DA) == [
    "2set_features.json",
    "FeatureSetEnrichment",
    "feature_x_metric_x_score.tsv",
    "feature_x_sample_x_number.tsv",
    "set_features.json",
    "target_x_sample_x_number.tsv",
]

# ---- #

const SE = joinpath(DA, "set_features.json")

@test BioLab.Error.@is GSEA._read_set(SE, Vector{String}(), 33, 36, 0)

# ---- #

const SE1_, FE11___ = GSEA._read_set(SE, unique(vcat(values(BioLab.Dict.read(SE))...)), 33, 36, 0)

@test length(SE1_) === 2

@test length(FE11___) === 2

# ---- #

const SE2_, FE12___ = GSEA._read_set(SE, ["SHH", "XIST"], 1, 5656, 0)

@test length(SE2_) === 2

@test length(FE12___) === 2

# ---- #

for (al, re) in (
    ("ks", GSEA.FeatureSetEnrichment.KS()),
    ("ksa", GSEA.FeatureSetEnrichment.KSa()),
    ("kli", GSEA.FeatureSetEnrichment.KLi()),
    ("kliom", GSEA.FeatureSetEnrichment.KLioM()),
    ("kliop", GSEA.FeatureSetEnrichment.KLioP()),
)

    GSEA._set_algorithm(al) == re

end

# ---- #

const DAO = joinpath(dirname(@__DIR__), "benchmark", "outer_loop")

const COT = joinpath(TE, "target_x_sample_x_number.tsv")

const COF = joinpath(TE, "feature_x_sample_x_number.tsv")

@test GSEA.convert_cls_gct(
    COT,
    COF,
    joinpath(DAO, "Coller_et_al_phen.cls"),
    joinpath(DAO, "Coller_et_al_gene_exp_preproc.gct"),
) === (COT, COF)

# ---- #

@test read(COF) == read(joinpath(DAO, "input", "feature_x_sample_x_number.tsv"))

@test read(COT) == read(joinpath(DAO, "input", "target_x_sample_x_number.tsv"))

# ---- #

const COS = joinpath(TE, "set_features.json")

@test GSEA.convert_gmt(
    COS,
    joinpath(DAO, "h.all.v2022.1.Hs.symbols.gmt"),
    joinpath(DAO, "c2.all.v2022.1.Hs.symbols.gmt"),
) === COS

# ---- #

@test read(COS) == read(joinpath(DAO, "input", "set_features.json"))

# ---- #

const TSF = joinpath(DA, "feature_x_sample_x_number.tsv")

@test BioLab.Error.@is GSEA.data_rank("", TSF, SE)

# ---- #

const OUD = GSEA.data_rank(BioLab.Path.remake_directory(joinpath(TE, "data_rank")), TSF, SE)

# ---- #

const SET_X_SAMPLE_X_ENRICHMENT =
    BioLab.DataFrame.read(joinpath(OUD, "set_x_sample_x_enrichment.tsv"))

@test size(SET_X_SAMPLE_X_ENRICHMENT) === (8, 10)

@test SET_X_SAMPLE_X_ENRICHMENT[!, "Set"] == [
    "HALLMARK_ESTROGEN_RESPONSE_LATE",
    "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
    "HALLMARK_ESTROGEN_RESPONSE_EARLY",
    "HALLMARK_KRAS_SIGNALING_DN",
    "HALLMARK_IL2_STAT5_SIGNALING",
    "HALLMARK_APICAL_JUNCTION",
    "HALLMARK_HYPOXIA",
    "HALLMARK_GLYCOLYSIS",
]

# ---- #

const TSM = joinpath(DA, "feature_x_metric_x_score.tsv")

@test BioLab.Error.@is GSEA.user_rank("", TSM, SE)

# ---- #

const OUU = GSEA.user_rank(
    BioLab.Path.remake_directory(joinpath(TE, "user_rank")),
    TSM,
    SE;
    number_of_sets_to_plot = 2,
    more_sets_to_plot = "HALLMARK_MYC_TARGETS_V1 HALLMARK_UV_RESPONSE_DN HALLMARK_UV_RESPONSE_UP ALIEN",
)

# ---- #

const SET_X_STATISTIC_X_NUMBER =
    BioLab.DataFrame.read(joinpath(OUU, "set_x_statistic_x_number.tsv"))

@test size(SET_X_STATISTIC_X_NUMBER) === (50, 5)

@test names(SET_X_STATISTIC_X_NUMBER) ==
      ["Set", "Enrichment", "Normalized Enrichment", "P-Value", "Adjusted P-Value"]

# ---- #

for (id, re) in (
    (1, ["HALLMARK_PANCREAS_BETA_CELLS", -0.35266, -1.36616, 0.0200837]),
    (2, ["HALLMARK_PROTEIN_SECRETION", -0.272096, -1.25207, 0.0686192]),
    (49, ["HALLMARK_MYC_TARGETS_V1", 0.603356, 2.73998, 0.000262812]),
    (50, ["HALLMARK_MYC_TARGETS_V2", 0.866579, 3.36557, 0.000262812]),
)

    @test SET_X_STATISTIC_X_NUMBER[id, 1] === re[1]

    @test isapprox(collect(SET_X_STATISTIC_X_NUMBER[id, 2:length(re)]), re[2:end]; atol = 1e-5)

end

@test length(BioLab.Path.read(OUU; ke_ = (r"html$",))) === 6

# ---- #

for (nu1_, nu2_, re) in (
    (fill(0, 3), fill(0, 3), 0.0),
    (fill(1, 3), fill(0.001, 3), 4.990009990009989),
    (collect(1:3), collect(10:10:30), -1.6363636363636365),
    (fill(0.001, 3), fill(1, 3), -4.990009990009989),
    (fill(0.001, 3), fill(10, 3), -4.999000099990002),
)

    @test GSEA._get_signal_to_noise_ratio(nu1_, nu2_) === re

    # 48.196 ns (0 allocations: 0 bytes)
    # 21.899 ns (0 allocations: 0 bytes)
    # 57.646 ns (0 allocations: 0 bytes)
    # 57.631 ns (0 allocations: 0 bytes)
    # 57.587 ns (0 allocations: 0 bytes)
    #@btime GSEA._get_signal_to_noise_ratio($nu1_, $nu2_)

end

# ---- #

const TST = joinpath(DA, "target_x_sample_x_number.tsv")

@test BioLab.Error.@is GSEA.metric_rank("", TST, TSF, SE)

# ---- #

const OUM =
    GSEA.metric_rank(BioLab.Path.remake_directory(joinpath(TE, "metric_rank")), TST, TSF, SE)

# ---- #

const FEATURE_X_METRIC_X_SCORE =
    BioLab.DataFrame.read(joinpath(OUM, "feature_x_metric_x_score.tsv"))

@test size(FEATURE_X_METRIC_X_SCORE) === (1000, 2)

@test names(FEATURE_X_METRIC_X_SCORE) == ["Feature", "signal-to-noise-ratio"]

@test isapprox(view(FEATURE_X_METRIC_X_SCORE, [1, 1000], 2), [1.7411, -1.83724]; atol = 1e-5)

# ---- #

const SET_X_STATISTIC_X_NUMBER =
    BioLab.DataFrame.read(joinpath(OUM, "set_x_statistic_x_number.tsv"))

@test size(SET_X_STATISTIC_X_NUMBER) === (8, 5)

@test length(BioLab.Path.read(OUM; ke_ = (r"html$",))) === 8

# ---- #

@test length(
    BioLab.Path.read(
        GSEA.metric_rank(
            BioLab.Path.remake_directory(joinpath(TE, "metric_rank_small")),
            TST,
            TSF,
            joinpath(DA, "2set_features.json");
            minimum_set_size = 3,
            maximum_set_size = 3,
            number_of_sets_to_plot = 2,
        );
        ke_ = (r"html$",),
    ),
) === 2
