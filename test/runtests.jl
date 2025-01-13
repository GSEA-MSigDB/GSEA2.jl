using GSEA

using Test: @test

# ----------------------------------------------------------------------------------------------- #

using Random: seed!

using Omics

# ---- #

const SC = -2.0

# 3.666 ns (0 allocations: 0 bytes)
# 1.750 ns (0 allocations: 0 bytes)
# 11.386 ns (0 allocations: 0 bytes)
# 11.386 ns (0 allocations: 0 bytes)
# 1.458 ns (0 allocations: 0 bytes)
# 3.958 ns (0 allocations: 0 bytes)
for ex in (-1, 0, 0.0, 0.1, 0.5, 1, 2)

    @test GSEA._absolute_exponentiate(SC, ex) === abs(SC)^ex

    @btime GSEA._absolute_exponentiate(SC, $ex)

    @btime abs(SC)^$ex

end

# ---- #

const SC_ = [-2, -1, -0.5, 0, 0, 0.5, 1, 2, 3.4]

const UF = lastindex(SC_)

# TODO: Benchmark Vector{Bool}
const IS_ = BitVector((1, 0, 1, 0, 1, 1, 0, 0, 1))

# ---- #

# 63.350 ns (0 allocations: 0 bytes)
# 63.350 ns (0 allocations: 0 bytes)
# 63.350 ns (0 allocations: 0 bytes)
# 9.510 ns (0 allocations: 0 bytes)
# 20.687 ns (0 allocations: 0 bytes)
for (ex, re) in (
    (-0.5, (UF, 0.0)),
    (0.1, (UF, 0.24581982412836917)),
    (0.5, (UF, 0.21402570288861142)),
    (1, (UF, 0.15625)),
    (2, (UF, 0.06226650062266501)),
)

    @test GSEA._get_1_normalizer(SC_, ex, IS_) === re

    @btime GSEA._get_1_normalizer(SC_, $ex, IS_)

end

# ---- #

const U2 = 0.25

# 59.700 ns (0 allocations: 0 bytes)
# 59.682 ns (0 allocations: 0 bytes)
# 59.681 ns (0 allocations: 0 bytes)
# 9.500 ns (0 allocations: 0 bytes)
# 21.314 ns (0 allocations: 0 bytes)
for (ex, re) in (
    (-0.5, (UF, U2, 0.0)),
    (0.1, (UF, U2, 0.24581982412836917)),
    (0.5, (UF, U2, 0.21402570288861142)),
    (1, (UF, U2, 0.15625)),
    (2, (UF, U2, 0.06226650062266501)),
)

    @test GSEA._get_0_1_normalizer(SC_, ex, IS_) === re

    @btime GSEA._get_0_1_normalizer(SC_, $ex, IS_)

end

# ---- #

# 97.238 ns (0 allocations: 0 bytes)
# 97.238 ns (0 allocations: 0 bytes)
# 97.281 ns (0 allocations: 0 bytes)
# 8.884 ns (0 allocations: 0 bytes)
# 28.098 ns (0 allocations: 0 bytes)
for (ex, re) in (
    (-0.5, (UF, 0.0, 0.0)),
    (0.1, (UF, 0.14006007078470165, 0.24581982412836917)),
    (0.5, (UF, 0.12366213677204271, 0.21402570288861142)),
    (1, (UF, 0.09615384615384615, 0.15625)),
    (2, (UF, 0.04533091568449683, 0.06226650062266501)),
)

    @test GSEA._get_all_1_normalizer(SC_, ex, IS_) === re

    @btime GSEA._get_all_1_normalizer(SC_, $ex, IS_)

end

# ---- #

# 1.458 ns (0 allocations: 0 bytes)
for (noa, no1, re) in ((0.5, 1 / 3, -1.0),)

    @test GSEA._get_0_normalizer(noa, no1) === re

    @btime GSEA._get_0_normalizer($noa, $no1)

end

# ---- #

const CA_ = ["K", "Q", "J", "X", "9", "8", "7", "6", "5", "4", "3", "2", "A"]

const CR_ = [6.0, 5, 4, 3, 2, 1, 0, -1, -2, -3, -4, -5, -6]

const CD_ = ["K", "A"]

const CC_ = map(in(Set(CD_)), CA_)

# ---- #

# 19.433 ns (0 allocations: 0 bytes)
# 17.869 ns (0 allocations: 0 bytes)
# 118.359 ns (0 allocations: 0 bytes)
# 126.633 ns (0 allocations: 0 bytes)
# 225.103 ns (0 allocations: 0 bytes)
# 225.069 ns (0 allocations: 0 bytes)
for (al, re) in zip(AL_, (-0.5, 0.0, 0.0, 0.0, 0.0, 0.0))

    @test isapprox(GSEA._enrich!(al, CR_, EX, CC_, nothing), re; atol=1e-15)

    @btime GSEA._enrich!($al, CR_, EX, CC_, nothing)

end

# ---- #

for al in AL_

    GSEA.plot("", al, CA_, CR_, CD_; title_text = al)

end

# ---- #

include("get_normalizer.jl")

# ---- #

const AL_ = GSEA.KS(), GSEA.KSa(), GSEA.KLi1(), GSEA.KLi(), GSEA.KLioM(), GSEA.KLioP()

const EX = 1

# ---- #

include("card.jl")

# ---- #

const DA = pkgdir(GSEA, "data")

# ---- #

const FE_, SO_ = eachcol(
    reverse!(
        Nucleus.DataFrame.read(
            joinpath(DA, "gene_x_statistic_x_number.tsv");
            select = [1, 2],
        ),
    ),
)

# ---- #

const FE1_ =
    Nucleus.GMT.read(joinpath(DA, "c2.all.v7.1.symbols.gmt"))["COLLER_MYC_TARGETS_UP"]

# ---- #

for al in AL_

    GSEA.plot("", al, FE_, SO_, FE1_; ex = EX, title_text = GSEA.make_string(al))

end

# ---- #

const IS_ = in(Set(FE1_)).(FE_)

# ---- #

# 45.208 μs (0 allocations: 0 bytes)
# 37.500 μs (0 allocations: 0 bytes)
# 164.833 μs (0 allocations: 0 bytes)
# 186.208 μs (0 allocations: 0 bytes)
# 325.833 μs (0 allocations: 0 bytes)
# 326.042 μs (0 allocations: 0 bytes)
for (al, re) in zip(
    AL_,
    (
        0.7651927829281453,
        0.41482514169516305,
        0.8524266036047564,
        0.7736480596525319,
        0.7750661968892066,
        0.772229922415844,
    ),
)

    @test isapprox(GSEA._enrich!(al, SO_, EX, IS_, nothing), re; atol = 0.000000000001)

    @btime GSEA._enrich!($al, SO_, EX, IS_, nothing)

end

# ---- #

const SE_FE1_ = Nucleus.GMT.read(joinpath(DA, "h.all.v7.1.symbols.gmt"))

# ---- #

const SE_ = collect(keys(SE_FE1_))

# ---- #

const FE1___ = collect(values(SE_FE1_))

# ---- #

# 3.022 ms (108 allocations: 934.22 KiB)
# 2.649 ms (108 allocations: 934.22 KiB)
# 9.001 ms (108 allocations: 934.22 KiB)
# 10.127 ms (108 allocations: 934.22 KiB)
# 17.202 ms (108 allocations: 934.22 KiB)
# 17.186 ms (108 allocations: 934.22 KiB)
for al in AL_

    @btime GSEA.enrich($al, FE_, SO_, FE1___; ex = EX)

end

# ---- #

const FE_X_SA_X_SC = hcat(SO_, SO_ * 10, fill(0.8, lastindex(FE_)))

# ---- #

# 9.548 ms (370 allocations: 5.51 MiB)
# 8.382 ms (370 allocations: 5.51 MiB)
# 27.547 ms (370 allocations: 5.51 MiB)
# 30.802 ms (370 allocations: 5.51 MiB)
# 52.180 ms (370 allocations: 5.51 MiB)
# 52.183 ms (370 allocations: 5.51 MiB)
for al in AL_

    @btime GSEA.enrich($al, FE_, FE_X_SA_X_SC, FE1___; ex = EX)

end

# ---- #

const AL = GSEA.KS()

# ---- #

const SE_X_SA_X_EN = GSEA.enrich(AL, FE_, FE_X_SA_X_SC, FE1___; ex = EX)

# ---- #

const TE = tempdir()

# ---- #

GSEA.plot(
    TE,
    AL,
    FE_,
    FE_X_SA_X_SC,
    FE1___,
    "Sample",
    SE_,
    ["Score", "Score x 10", "Constant"],
    SE_X_SA_X_EN;
    ex = EX,
)

# ---- #

const JS = joinpath(DA, "set_features.json")

# ---- #

for (fe_, mi, ma, fr) in (
    (String[], 33, 36, 0),
    (unique(vcat(values(Nucleus.Dict.read(JS))...)), 33, 36, 0),
    (["SHH", "XIST"], 1, 5656, 0),
)

    if isempty(fe_)

        @test Nucleus.Error.@is GSEA._read_set(JS, fe_, mi, ma, fr)

    else

        se_, fe1___ = GSEA._read_set(JS, fe_, mi, ma, fr)

        @test lastindex(se_) === lastindex(fe1___) === 2

    end

end

# ---- #

for (al, re) in zip(("ks", "ksa", "kli", "kliom", "kliop"), AL_)

    GSEA._set_algorithm(al) == re

end

# ---- #

const TSF = joinpath(DA, "feature_x_sample_x_number.tsv")

# ---- #

const OUD = mkdir(joinpath(TE, "data_rank"))

# ---- #

GSEA.data_rank(OUD, TSF, JS)

# ---- #

const SET_X_SAMPLE_X_ENRICHMENT =
    Nucleus.DataFrame.read(joinpath(OUD, "set_x_sample_x_enrichment.tsv"))

# ---- #

@test size(SET_X_SAMPLE_X_ENRICHMENT) === (8, 10)

# ---- #

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

# 737.209 μs (1001 allocations: 1.47 MiB)
# 803.333 μs (1001 allocations: 1.47 MiB)
for al in (AL_[1], AL_[end])

    seed!(20231103)

    en_ = randn(100)

    se_x_id_x_ra = randn(100, 1000)

    GSEA._normalize_enrichment!(al, en_, se_x_id_x_ra)

    @btime GSEA._normalize_enrichment!($al, $en_, $se_x_id_x_ra)

end

# ---- #

function test_statistics(set_x_statistic_x_numberu, n_ro)

    @test size(set_x_statistic_x_numberu, 1) === n_ro

    @test names(set_x_statistic_x_numberu) ==
          ["Set", "Enrichment", "Normalized Enrichment", "P-Value", "Adjusted P-Value"]

end

# ---- #

function test_html(ou, uh)

    @test lastindex(Nucleus.Path.read(ou; ke_ = (r"html$",))) === uh

end

# ---- #

const OUU = mkdir(joinpath(TE, "user_rank"))

# ---- #

GSEA.user_rank(
    OUU,
    joinpath(DA, "feature_x_metric_x_score.tsv"),
    JS;
    number_of_sets_to_plot = 2,
    more_sets_to_plot = "HALLMARK_MYC_TARGETS_V1 HALLMARK_UV_RESPONSE_DN HALLMARK_UV_RESPONSE_UP ALIEN",
)

# ---- #

const SET_X_STATISTIC_X_NUMBERU =
    Nucleus.DataFrame.read(joinpath(OUU, "set_x_statistic_x_number.tsv"))

# ---- #

test_statistics(SET_X_STATISTIC_X_NUMBERU, 50)

# ---- #

for (id, re1, re2) in (
    (1, "HALLMARK_PANCREAS_BETA_CELLS", [-0.35266, -1.36616, 0.0200837]),
    (2, "HALLMARK_PROTEIN_SECRETION", [-0.272096, -1.25207, 0.0686192]),
    (49, "HALLMARK_MYC_TARGETS_V1", [0.603356, 2.73998, 0.000262812]),
    (50, "HALLMARK_MYC_TARGETS_V2", [0.866579, 3.36557, 0.000262812]),
)

    @test SET_X_STATISTIC_X_NUMBERU[id, 1] === re1

    @test isapprox(collect(SET_X_STATISTIC_X_NUMBERU[id, 2:4]), re2; atol = 0.00001)

end

# ---- #

test_html(OUU, 6)

# ---- #

const ZE_ = fill(0, 3)

# ---- #

@test isnan(GSEA._get_signal_to_noise_ratio(ZE_, ZE_))

# ---- #

# 37.088 ns (0 allocations: 0 bytes)
# 21.439 ns (0 allocations: 0 bytes)
# 38.807 ns (0 allocations: 0 bytes)
# 39.396 ns (0 allocations: 0 bytes)
for (nu1_, nu2_, re) in (
    (fill(1, 3), fill(0.001, 3), 4.990009990009989),
    (collect(1:3), collect(10:10:30), -1.6363636363636365),
    (fill(0.001, 3), fill(1, 3), -4.990009990009989),
    (fill(0.001, 3), fill(10, 3), -4.999000099990002),
)

    @test GSEA._get_signal_to_noise_ratio(nu1_, nu2_) === re

    @btime GSEA._get_signal_to_noise_ratio($nu1_, $nu2_)

end

# ---- #

const TST = joinpath(DA, "target_x_sample_x_number.tsv")

# ---- #

const OUM = mkdir(joinpath(TE, "metric_rank"))

# ---- #

GSEA.metric_rank(OUM, TST, TSF, JS)

# ---- #

const FEATURE_X_METRIC_X_SCORE =
    Nucleus.DataFrame.read(joinpath(OUM, "feature_x_metric_x_score.tsv"))

# ---- #

@test size(FEATURE_X_METRIC_X_SCORE) === (1000, 2)

# ---- #

@test names(FEATURE_X_METRIC_X_SCORE) == ["Feature", "signal-to-noise-ratio"]

# ---- #

@test isapprox(
    view(FEATURE_X_METRIC_X_SCORE, [1, 1000], 2),
    [1.83724, -1.7411];
    atol = 0.00001,
)

# ---- #

const SET_X_STATISTIC_X_NUMBERM =
    Nucleus.DataFrame.read(joinpath(OUM, "set_x_statistic_x_number.tsv"))

# ---- #

test_statistics(SET_X_STATISTIC_X_NUMBERM, 8)

# ---- #

test_html(OUM, 8)

# ---- #

const OUMS = mkdir(joinpath(TE, "metric_rank_small"))

# ---- #

GSEA.metric_rank(
    OUMS,
    TST,
    TSF,
    joinpath(DA, "2set_features.json");
    minimum_set_size = 3,
    maximum_set_size = 3,
    number_of_sets_to_plot = 2,
)

# ---- #

test_html(OUMS, 2)
