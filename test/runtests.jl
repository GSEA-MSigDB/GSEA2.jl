using GSEA

using Test: @test

# ----------------------------------------------------------------------------------------------- #

using Omics

using CLSGCTGMT

# ---- #

const SC_ = [-2, -1, -0.5, 0, 0, 0.5, 1, 2, 3.4]

const BO_ = [true, false, true, false, true, true, false, false, true]

# ---- #

const N0 = -0.25

# 59.682 ns (0 allocations: 0 bytes)
# 59.681 ns (0 allocations: 0 bytes)
# 9.500 ns (0 allocations: 0 bytes)
# 21.314 ns (0 allocations: 0 bytes)
#
# 69.928 ns (0 allocations: 0 bytes)
# 69.928 ns (0 allocations: 0 bytes)
# 8.291 ns (0 allocations: 0 bytes)
# 19.831 ns (0 allocations: 0 bytes)
for (ex, re) in (
    (0.1, (N0, 0.24581982412836917)),
    (0.5, (N0, 0.21402570288861142)),
    (1, (N0, 0.15625)),
    (2, (N0, 0.06226650062266501)),
)

    @test GSEA._get_0_1(SC_, ex, BO_) === re

    #@btime GSEA._get_0_1(SC_, $ex, BO_)

end

# ---- #

# 97.238 ns (0 allocations: 0 bytes)
# 97.281 ns (0 allocations: 0 bytes)
# 8.884 ns (0 allocations: 0 bytes)
# 28.098 ns (0 allocations: 0 bytes)
#
# 116.858 ns (0 allocations: 0 bytes)
# 116.903 ns (0 allocations: 0 bytes)
# 8.291 ns (0 allocations: 0 bytes)
# 24.992 ns (0 allocations: 0 bytes)
for (ex, re) in (
    (0.1, (0.14006007078470165, 0.24581982412836917)),
    (0.5, (0.12366213677204271, 0.21402570288861142)),
    (1, (0.09615384615384615, 0.15625)),
    (2, (0.04533091568449683, 0.06226650062266501)),
)

    @test GSEA._get_all_1(SC_, ex, BO_) === re

    #@btime GSEA._get_all_1(SC_, $ex, BO_)

end

# ---- #

# 1.458 ns (0 allocations: 0 bytes)
#
# 2.083 ns (0 allocations: 0 bytes)
for (sa, s1, re) in ((0.5, 1 / 3, -1.0),)

    @test GSEA._get_0(sa, s1) === re

    #@btime GSEA._get_0($sa, $s1)

end

# ---- #

# 63.350 ns (0 allocations: 0 bytes)
# 63.350 ns (0 allocations: 0 bytes)
# 9.510 ns (0 allocations: 0 bytes)
# 20.687 ns (0 allocations: 0 bytes)
#
# 79.287 ns (0 allocations: 0 bytes)
# 79.295 ns (0 allocations: 0 bytes)
# 10.750 ns (0 allocations: 0 bytes)
# 17.242 ns (0 allocations: 0 bytes)
for (ex, re) in (
    (0.1, (0.24581982412836917)),
    (0.5, (0.21402570288861142)),
    (1, (0.15625)),
    (2, (0.06226650062266501)),
)

    @test GSEA._get_1(SC_, ex, BO_) === re

    #@btime GSEA._get_1(SC_, $ex, BO_)

end

# ---- #

const AL_ = GSEA.KS(), GSEA.KSa(), GSEA.KLioM(), GSEA.KLioP(), GSEA.KLi(), GSEA.KLi1()

# ---- #

const DA = pkgdir(GSEA, "data")

# ---- #

const FE_, SO_ = eachcol(
    reverse!(
        Omics.Table.rea(joinpath(DA, "gene_x_statistic_x_number.tsv"); select = [1, 2]),
    ),
)

# ---- #

# 19.433 ns (0 allocations: 0 bytes)
# 17.869 ns (0 allocations: 0 bytes)
# 225.103 ns (0 allocations: 0 bytes)
# 225.069 ns (0 allocations: 0 bytes)
# 126.633 ns (0 allocations: 0 bytes)
# 118.359 ns (0 allocations: 0 bytes)
# 45.208 μs (0 allocations: 0 bytes)
# 37.500 μs (0 allocations: 0 bytes)
# 325.833 μs (0 allocations: 0 bytes)
# 326.042 μs (0 allocations: 0 bytes)
# 186.208 μs (0 allocations: 0 bytes)
# 164.833 μs (0 allocations: 0 bytes)
#
# 239.904 ns (6 allocations: 400 bytes)
# 17.285 ns (0 allocations: 0 bytes)
# 16.658 ns (0 allocations: 0 bytes)
# 285.187 ns (0 allocations: 0 bytes)
# 285.187 ns (0 allocations: 0 bytes)
# 155.514 ns (0 allocations: 0 bytes)
# 143.777 ns (0 allocations: 0 bytes)
# 460.000 μs (7 allocations: 20.42 KiB)
# 43.375 μs (0 allocations: 0 bytes)
# 37.541 μs (0 allocations: 0 bytes)
# 413.083 μs (0 allocations: 0 bytes)
# 413.166 μs (0 allocations: 0 bytes)
# 243.416 μs (0 allocations: 0 bytes)
# 210.291 μs (0 allocations: 0 bytes)
for (fe_, sc_, me_, re_) in (
    (
        ["K", "Q", "J", "X", "9", "8", "7", "6", "5", "4", "3", "2", "A"],
        [6.0, 5, 4, 3, 2, 1, 0, -1, -2, -3, -4, -5, -6],
        ["K", "A"],
        (-0.5, 0.0, 0.0, 0.0, 0.0, 0.0),
    ),
    (
        FE_,
        SO_,
        CLSGCTGMT.read_gmt(joinpath(DA, "c2.all.v7.1.symbols.gmt"))["COLLER_MYC_TARGETS_UP"],
        (
            0.7651927829281453,
            0.41482514169516305,
            1.1181841586127337,
            1.1140922794954267,
            1.1161382190540838,
            1.2297916337424049,
        ),
    ),
)

    bo_ = GSEA._is_in(fe_, me_)

    @test typeof(bo_) === Vector{Bool}

    #@btime GSEA._is_in($fe_, $me_)

    # TODO
    GSEA._is_in!

    for (al, re) in zip(AL_, re_)

        ex = 1.0

        @test isapprox(GSEA._enrich!(al, sc_, ex, bo_, nothing), re; atol = 1e-11)

        #@btime GSEA._enrich!($al, $sc_, $ex, $bo_, nothing)

        GSEA.plot(
            "",
            al,
            fe_,
            sc_,
            me_;
            la = Dict("title" => Dict("text" => string(al)[6:(end - 2)])),
        )

    end

end

# ---- #

const SE_ME_ = CLSGCTGMT.read_gmt(joinpath(DA, "h.all.v7.1.symbols.gmt"))

const SE_ = collect(keys(SE_ME_))

const ME___ = collect(values(SE_ME_))

const UE = lastindex(SE_)

# ---- #

# 3.022 ms (108 allocations: 934.22 KiB)
# 2.649 ms (108 allocations: 934.22 KiB)
# 17.202 ms (108 allocations: 934.22 KiB)
# 17.186 ms (108 allocations: 934.22 KiB)
# 10.127 ms (108 allocations: 934.22 KiB)
# 9.001 ms (108 allocations: 934.22 KiB)
#
# 2.945 ms (13 allocations: 803.23 KiB)
# 2.683 ms (13 allocations: 803.23 KiB)
# 21.294 ms (13 allocations: 803.23 KiB)
# 21.374 ms (13 allocations: 803.23 KiB)
# 13.060 ms (13 allocations: 803.23 KiB)
# 11.431 ms (13 allocations: 803.23 KiB)
for al in AL_

    @test lastindex(GSEA.enrich(al, FE_, SO_, ME___)) === UE

    #@btime GSEA.enrich($al, FE_, SO_, ME___)

end

# ---- #

const SC = hcat(SO_, SO_ * 10.0, fill(0.8, lastindex(FE_)))

const RE = UE, size(SC, 2)

# 9.548 ms (370 allocations: 5.51 MiB)
# 8.382 ms (370 allocations: 5.51 MiB)
# 52.180 ms (370 allocations: 5.51 MiB)
# 52.183 ms (370 allocations: 5.51 MiB)
# 30.802 ms (370 allocations: 5.51 MiB)
# 27.547 ms (370 allocations: 5.51 MiB)
#
# 9.463 ms (99 allocations: 6.03 MiB)
# 8.600 ms (99 allocations: 6.03 MiB)
# 64.487 ms (99 allocations: 6.03 MiB)
# 64.616 ms (99 allocations: 6.03 MiB)
# 39.478 ms (99 allocations: 6.03 MiB)
# 34.835 ms (99 allocations: 6.03 MiB)
for al in AL_

    @test size(GSEA.enrich(al, FE_, SC, ME___)) === RE

    #@btime GSEA.enrich($al, FE_, SC, ME___)

end

# ---- #

const JS = joinpath(DA, "set_features.json")

const OU = joinpath(homedir(), "Downloads")

# ---- #

const AL = GSEA.KS()

GSEA.plot(
    tempdir(),
    FE_,
    SC,
    AL,
    SE_,
    ME___,
    "Sample",
    ["Score", "Score x 10", "Constant"],
    GSEA.enrich(AL, FE_, SC, ME___),
)

# ---- #

for (fe_, mi, ma, fr) in (
    (String[], 33, 36, 0),
    (unique(vcat(values(Omics.Dic.rea(JS))...)), 33, 36, 0),
    (["SHH", "XIST"], 1, 5656, 0),
)

    if isempty(fe_)

        @test Omics.Error.@is GSEA._read_set(JS, fe_, mi, ma, fr)

    else

        se_, me___ = GSEA._read_set(JS, fe_, mi, ma, fr)

        @test lastindex(se_) === lastindex(me___) === 2

    end

end

# ---- #

for (al, re) in zip(("ks", "ksa", "kli", "kliom", "kliop"), AL_)

    GSEA._set_algorithm(al) == re

end

# ---- #

const TF = joinpath(DA, "feature_x_sample_x_number.tsv")

# ---- #

const OD = mkpath(joinpath(OU, "data_rank"))

GSEA.data_rank(OD, TF, JS)

const TD = Omics.Table.rea(joinpath(OD, "TD.tsv"))

@test size(TD) === (8, 10)

@test TD[!, "Set"] == [
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

    ra = randn(100, 1000)

    GSEA._normalize_enrichment!(al, en_, ra)

    @btime GSEA._normalize_enrichment!($al, $en_, $ra)

end

# ---- #

function test_statistics(ta, ur)

    @test size(ta, 1) === ur

    @test names(ta) == ["Set", "Enrichment", "Normalized Enrichment", "P-Value", "Q-Value"]

end

function test_html(ou, uh)

    @test lastindex(readdir(ou; ke_ = (r"html$",))) === uh

end

# ---- #

const OS = mkpath(joinpath(OU, "user_rank"))

GSEA.user_rank(
    OS,
    joinpath(DA, "feature_x_metric_x_score.tsv"),
    JS;
    number_of_sets_to_plot = 2,
    more_sets_to_plot = "HALLMARK_MYC_TARGETS_V1 HALLMARK_UV_RESPONSE_DN HALLMARK_UV_RESPONSE_UP ALIEN",
)

const SU = Omics.Table.read(joinpath(OS, "set_x_statistic_x_number.tsv"))

test_statistics(SU, 50)

for (id, r1, r2) in (
    (1, "HALLMARK_PANCREAS_BETA_CELLS", [-0.35266, -1.36616, 0.0200837]),
    (2, "HALLMARK_PROTEIN_SECRETION", [-0.272096, -1.25207, 0.0686192]),
    (49, "HALLMARK_MYC_TARGETS_V1", [0.603356, 2.73998, 0.000262812]),
    (50, "HALLMARK_MYC_TARGETS_V2", [0.866579, 3.36557, 0.000262812]),
)

    @test SU[id, 1] === r1

    @test isapprox(collect(SU[id, 2:4]), r2; atol = 0.00001)

end

test_html(OS, 6)

# ---- #

const ZE_ = zeros(Int, 3)

@test isnan(GSEA._get_signal_to_noise_ratio(ZE_, ZE_))

# 37.088 ns (0 allocations: 0 bytes)
# 21.439 ns (0 allocations: 0 bytes)
# 38.807 ns (0 allocations: 0 bytes)
# 39.396 ns (0 allocations: 0 bytes)
for (n1_, n2_, re) in (
    (fill(1, 3), fill(0.001, 3), 4.990009990009989),
    (collect(1:3), collect(10:10:30), -1.6363636363636365),
    (fill(0.001, 3), fill(1, 3), -4.990009990009989),
    (fill(0.001, 3), fill(10, 3), -4.999000099990002),
)

    @test GSEA._get_signal_to_noise_ratio(n1_, n2_) === re

    @btime GSEA._get_signal_to_noise_ratio($n1_, $n2_)

end

# ---- #

const TT = joinpath(DA, "target_x_sample_x_number.tsv")

# ---- #

const OM = mkpath(joinpath(OU, "metric_rank"))

GSEA.metric_rank(OM, TT, TF, JS)

const TM = Omics.Table.read(joinpath(OM, "feature_x_metric_x_score.tsv"))

@test size(TM) === (1000, 2)

@test names(TM) == ["Feature", "signal-to-noise-ratio"]

@test isapprox(view(TM, [1, 1000], 2), [1.83724, -1.7411]; atol = 0.00001)

test_statistics(Omics.Table.read(joinpath(OM, "set_x_statistic_x_number.tsv")), 8)

test_html(OM, 8)

# ---- #

const TR = mkpath(joinpath(OU, "metric_rank_small"))

GSEA.metric_rank(
    TR,
    TT,
    TF,
    joinpath(DA, "2set_features.json");
    minimum_set_size = 3,
    maximum_set_size = 3,
    number_of_sets_to_plot = 2,
)

test_html(TR, 2)
