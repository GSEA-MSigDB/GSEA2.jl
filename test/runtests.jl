using Test: @test

using GSEA

# ----------------------------------------------------------------------------------------------- #

using Nucleus

# ---- #

const DA = joinpath(dirname(@__DIR__), "data")

# ---- #

@test Nucleus.Path.read(DA) == [
    "2set_features.json",
    "c2.all.v7.1.symbols.gmt",
    "feature_x_metric_x_score.tsv",
    "feature_x_sample_x_number.tsv",
    "gene_x_statistic_x_number.tsv",
    "h.all.v7.1.symbols.gmt",
    "set_features.json",
    "target_x_sample_x_number.tsv",
]

# ---- #

const AL_ = (GSEA.KS(), GSEA.KSa(), GSEA.KLi1(), GSEA.KLi(), GSEA.KLioM(), GSEA.KLioP())

# ---- #

for (al, re) in zip(AL_, ("KS", "KSa", "KLi1", "KLi", "KLioM", "KLioP"))

    @test GSEA.make_string(al) === re

    # 408.750 ns (9 allocations: 360 bytes)
    # 404.165 ns (9 allocations: 360 bytes)
    # 414.784 ns (9 allocations: 360 bytes)
    # 406.875 ns (9 allocations: 360 bytes)
    # 417.296 ns (9 allocations: 360 bytes)
    # 418.131 ns (9 allocations: 360 bytes)
    #@btime GSEA.make_string($al)

end

# ---- #

const SC = -2.0

# ---- #

for ex in (-1, 0, 1, 2, 3, 4, 0.1, 0.5)

    @test GSEA._absolute_exponentiate(SC, ex) === abs(SC)^ex

    # 3.666 ns (0 allocations: 0 bytes)
    # 1.750 ns (0 allocations: 0 bytes)
    # 1.458 ns (0 allocations: 0 bytes)
    # 3.958 ns (0 allocations: 0 bytes)
    # 2.375 ns (0 allocations: 0 bytes)
    # 4.875 ns (0 allocations: 0 bytes)
    # 11.386 ns (0 allocations: 0 bytes)
    # 11.386 ns (0 allocations: 0 bytes)
    #@btime GSEA._absolute_exponentiate(SC, $ex)

end

# ---- #

const SCS_ = [-2, -1, -0.5, 0, 0, 0.5, 1, 2, 3.4]

# ---- #

const N = lastindex(SCS_)

# ---- #

const ISS_ = BitVector((1, 0, 1, 0, 1, 1, 0, 0, 1))

# ---- #

for (ex, re) in (
    (-0.5, (N, 0.0)),
    (1, (N, 0.15625)),
    (2, (N, 0.06226650062266501)),
    (0.1, (N, 0.24581982412836917)),
    (0.5, (N, 0.21402570288861142)),
)

    @test GSEA._get_1_normalizer(SCS_, ex, ISS_) === re

    # 63.350 ns (0 allocations: 0 bytes)
    # 9.510 ns (0 allocations: 0 bytes)
    # 20.687 ns (0 allocations: 0 bytes)
    # 63.350 ns (0 allocations: 0 bytes)
    # 63.350 ns (0 allocations: 0 bytes)
    #@btime GSEA._get_1_normalizer(SCS_, $ex, ISS_)

end

# ---- #

const NOA = 0.25

# ---- #

for (ex, re) in (
    (-0.5, (N, NOA, 0.0)),
    (1, (N, NOA, 0.15625)),
    (2, (N, NOA, 0.06226650062266501)),
    (0.1, (N, NOA, 0.24581982412836917)),
    (0.5, (N, NOA, 0.21402570288861142)),
)

    @test GSEA._get_0_1_normalizer(SCS_, ex, ISS_) === re

    # 59.700 ns (0 allocations: 0 bytes)
    # 9.500 ns (0 allocations: 0 bytes)
    # 21.314 ns (0 allocations: 0 bytes)
    # 59.682 ns (0 allocations: 0 bytes)
    # 59.681 ns (0 allocations: 0 bytes)
    #@btime GSEA._get_0_1_normalizer(SCS_, $ex, ISS_)

end

# ---- #

for (ex, re) in (
    (-0.5, (N, 0.0, 0.0)),
    (1, (N, 0.09615384615384615, 0.15625)),
    (2, (N, 0.04533091568449683, 0.06226650062266501)),
    (0.1, (N, 0.14006007078470165, 0.24581982412836917)),
    (0.5, (N, 0.12366213677204271, 0.21402570288861142)),
)

    @test GSEA._get_all_1_normalizer(SCS_, ex, ISS_) === re

    # 97.238 ns (0 allocations: 0 bytes)
    # 8.884 ns (0 allocations: 0 bytes)
    # 28.098 ns (0 allocations: 0 bytes)
    # 97.238 ns (0 allocations: 0 bytes)
    # 97.281 ns (0 allocations: 0 bytes)
    #@btime GSEA._get_all_1_normalizer(SCS_, $ex, ISS_)

end

# ---- #

for (noa, no1, re) in ((0.5, 1 / 3, -1.0),)

    @test GSEA._get_0_normalizer(noa, no1) === re

    # 1.458 ns (0 allocations: 0 bytes)
    #@btime GSEA._get_0_normalizer($noa, $no1)

end

# ---- #

const FEC_ = ["K", "Q", "J", "X", "9", "8", "7", "6", "5", "4", "3", "2", "A"]

# ---- #

const SCC_ = [6.0, 5, 4, 3, 2, 1, 0, -1, -2, -3, -4, -5, -6]

# ---- #

const FE1C_ = ["K", "A"]

# ---- #

const EX = 1

# ---- #

for al in AL_

    GSEA.plot("", al, FEC_, SCC_, FE1C_; ex = EX, title_text = GSEA.make_string(al))

end

# ---- #

const ISC_ = in(Set(FE1C_)).(FEC_)

# ---- #

for (al, re) in zip(AL_, (-0.5, 0, 0, 0, 0, 0))

    @test isapprox(GSEA._enrich!(al, SCC_, EX, ISC_, nothing), re; atol = 0.000000000000001)

    # 19.433 ns (0 allocations: 0 bytes)
    # 17.869 ns (0 allocations: 0 bytes)
    # 118.359 ns (0 allocations: 0 bytes)
    # 126.633 ns (0 allocations: 0 bytes)
    # 225.103 ns (0 allocations: 0 bytes)
    # 225.069 ns (0 allocations: 0 bytes)
    #@btime GSEA._enrich!($al, SCC_, EX, ISC_, nothing)

end

# ---- #

const FE_, SC_ = eachcol(
    reverse!(
        Nucleus.DataFrame.read(joinpath(DA, "gene_x_statistic_x_number.tsv"); select = [1, 2]),
    ),
)

# ---- #

const FE1_ = Nucleus.GMT.read(joinpath(DA, "c2.all.v7.1.symbols.gmt"))["COLLER_MYC_TARGETS_UP"]

# ---- #

for al in AL_

    GSEA.plot("", al, FE_, SC_, FE1_; ex = EX, title_text = GSEA.make_string(al))

end

# ---- #

const IS_ = in(Set(FE1_)).(FE_)

# ---- #

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

    @test isapprox(GSEA._enrich!(al, SC_, EX, IS_, nothing), re; atol = 0.000000000001)

    # 45.208 μs (0 allocations: 0 bytes)
    # 37.500 μs (0 allocations: 0 bytes)
    # 164.833 μs (0 allocations: 0 bytes)
    # 186.208 μs (0 allocations: 0 bytes)
    # 325.833 μs (0 allocations: 0 bytes)
    # 326.042 μs (0 allocations: 0 bytes)
    #@btime GSEA._enrich!($al, SC_, EX, IS_, nothing)

end

# ---- #

const SE_FE1_ = Nucleus.GMT.read(joinpath(DA, "h.all.v7.1.symbols.gmt"))

# ---- #

const SE_ = collect(keys(SE_FE1_))

# ---- #

const FE1___ = collect(values(SE_FE1_))

# ---- #

for al in AL_

    # 3.022 ms (108 allocations: 934.22 KiB)
    # 2.649 ms (108 allocations: 934.22 KiB)
    # 9.001 ms (108 allocations: 934.22 KiB)
    # 10.139 ms (108 allocations: 934.22 KiB)
    # 17.202 ms (108 allocations: 934.22 KiB)
    # 17.186 ms (108 allocations: 934.22 KiB)
    #@btime GSEA.enrich($al, FE_, SC_, FE1___; ex = EX)

end

# ---- #

const FE_X_SA_X_SC = hcat(SC_, SC_ * 10, fill(0.8, lastindex(FE_)))

# ---- #

for al in AL_

    # 9.572 ms (370 allocations: 5.51 MiB)
    # 8.382 ms (370 allocations: 5.51 MiB)
    # 27.547 ms (370 allocations: 5.51 MiB)
    # 30.802 ms (370 allocations: 5.51 MiB)
    # 52.180 ms (370 allocations: 5.51 MiB)
    # 52.183 ms (370 allocations: 5.51 MiB)
    #@btime GSEA.enrich($al, FE_, FE_X_SA_X_SC, FE1___; ex = EX)

end

# ---- #

const AL = GSEA.KS()

# ---- #

const SE_X_SA_X_EN = GSEA.enrich(AL, FE_, FE_X_SA_X_SC, FE1___; ex = EX)

# ---- #

GSEA.plot(
    Nucleus.TE,
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

const TE = joinpath(tempdir(), "GSEA")

# ---- #

Nucleus.Path.remake_directory(TE)

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

const OUD = joinpath(TE, "data_rank")

# ---- #

Nucleus.Path.remake_directory(OUD)

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

function test_statistics(set_x_statistic_x_numberu, n_ro)

    @test size(set_x_statistic_x_numberu, 1) === n_ro

    @test names(set_x_statistic_x_numberu) ==
          ["Set", "Enrichment", "Normalized Enrichment", "P-Value", "Adjusted P-Value"]

end

# ---- #

function test_html(ou, n)

    @test lastindex(Nucleus.Path.read(ou; ke_ = (r"html$",))) === n

end

# ---- #

const OUU = joinpath(TE, "user_rank")

# ---- #

Nucleus.Path.remake_directory(OUU)

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

for (nu1_, nu2_, re) in (
    (fill(1, 3), fill(0.001, 3), 4.990009990009989),
    (collect(1:3), collect(10:10:30), -1.6363636363636365),
    (fill(0.001, 3), fill(1, 3), -4.990009990009989),
    (fill(0.001, 3), fill(10, 3), -4.999000099990002),
)

    @test GSEA._get_signal_to_noise_ratio(nu1_, nu2_) === re

    # 55.752 ns (0 allocations: 0 bytes)
    # 43.897 ns (0 allocations: 0 bytes)
    # 57.546 ns (0 allocations: 0 bytes)
    # 56.361 ns (0 allocations: 0 bytes)
    #@btime GSEA._get_signal_to_noise_ratio($nu1_, $nu2_)

end

# ---- #

const TST = joinpath(DA, "target_x_sample_x_number.tsv")

# ---- #

const OUM = joinpath(TE, "metric_rank")

# ---- #

Nucleus.Path.remake_directory(OUM)

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

@test isapprox(view(FEATURE_X_METRIC_X_SCORE, [1, 1000], 2), [1.83724, -1.7411]; atol = 0.00001)

# ---- #

const SET_X_STATISTIC_X_NUMBERM =
    Nucleus.DataFrame.read(joinpath(OUM, "set_x_statistic_x_number.tsv"))

# ---- #

test_statistics(SET_X_STATISTIC_X_NUMBERM, 8)

# ---- #

test_html(OUM, 8)

# ---- #

const OUMS = joinpath(TE, "metric_rank_small")

# ---- #

Nucleus.Path.remake_directory(OUMS)

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
