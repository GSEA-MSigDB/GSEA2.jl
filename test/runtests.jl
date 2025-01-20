using GSEA

using Test: @test

# ----------------------------------------------------------------------------------------------- #

using Random: seed!

using Omics

using CLSGCTGMT

# ---- #

# 16.073 ns (0 allocations: 0 bytes)
# 16.533 ns (0 allocations: 0 bytes)
# 16.367 ns (0 allocations: 0 bytes)
# 16.366 ns (0 allocations: 0 bytes)
for (n1_, n2_, re) in (
    (ones(Int, 3), fill(0.001, 3), 4.990009990009989),
    (collect(1:3), collect(10:10:30), -1.6363636363636365),
    (fill(0.001, 3), ones(Int, 3), -4.990009990009989),
    (fill(0.001, 3), fill(10, 3), -4.999000099990002),
)

    @test GSEA._get_signal_to_noise_ratio(n1_, n2_) === re

    #@btime GSEA._get_signal_to_noise_ratio($n1_, $n2_)

end

# ---- #

const AL_ = GSEA.KS(), GSEA.KSa(), GSEA.KLioM(), GSEA.KLioP(), GSEA.KLi(), GSEA.KLi1()

function strin(al)

    string(al)[6:(end - 2)]

end

# ---- #

const SC_ = [-2, -1, -0.5, 0, 0, 0.5, 1, 2, 3.4]

const BO_ = [true, false, true, false, true, true, false, false, true]

const N0 = -0.25

# 59.682 ns (0 allocations: 0 bytes)
# 59.681 ns (0 allocations: 0 bytes)
# 9.500 ns (0 allocations: 0 bytes)
# 21.314 ns (0 allocations: 0 bytes)
# 97.238 ns (0 allocations: 0 bytes)
# 97.281 ns (0 allocations: 0 bytes)
# 8.884 ns (0 allocations: 0 bytes)
# 28.098 ns (0 allocations: 0 bytes)
# 63.350 ns (0 allocations: 0 bytes)
# 63.350 ns (0 allocations: 0 bytes)
# 9.510 ns (0 allocations: 0 bytes)
# 20.687 ns (0 allocations: 0 bytes)
#
# 69.928 ns (0 allocations: 0 bytes)
# 69.928 ns (0 allocations: 0 bytes)
# 8.291 ns (0 allocations: 0 bytes)
# 19.747 ns (0 allocations: 0 bytes)
# 116.858 ns (0 allocations: 0 bytes)
# 116.857 ns (0 allocations: 0 bytes)
# 8.291 ns (0 allocations: 0 bytes)
# 25.016 ns (0 allocations: 0 bytes)
# 79.295 ns (0 allocations: 0 bytes)
# 79.287 ns (0 allocations: 0 bytes)
# 10.761 ns (0 allocations: 0 bytes)
# 17.242 ns (0 allocations: 0 bytes)
for (al, re_) in zip(
    AL_[[1, 3, 6]],
    (
        (
            (N0, 0.24581982412836917),
            (N0, 0.21402570288861142),
            (N0, 0.15625),
            (N0, 0.06226650062266501),
        ),
        (
            (0.14006007078470165, 0.24581982412836917),
            (0.12366213677204271, 0.21402570288861142),
            (0.09615384615384615, 0.15625),
            (0.04533091568449683, 0.06226650062266501),
        ),
        (0.24581982412836917, 0.21402570288861142, 0.15625, 0.06226650062266501),
    ),
)

    for (ex, re) in zip((0.1, 0.5, 1, 2), re_)

        @test GSEA._get_delta(al, SC_, ex, BO_) === re

        #@btime GSEA._get_delta($al, SC_, $ex, BO_)

    end

end

# ---- #

# 1.458 ns (0 allocations: 0 bytes)
#
# 2.083 ns (0 allocations: 0 bytes)
for (sa, s1, re) in ((0.5, 1 / 3, -1.0),)

    @test GSEA._get_delta(sa, s1) === re

    #@btime GSEA._get_delta($sa, $s1)

end

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
# 233.059 ns (6 allocations: 400 bytes)
# 17.285 ns (0 allocations: 0 bytes)
# 16.658 ns (0 allocations: 0 bytes)
# 284.452 ns (0 allocations: 0 bytes)
# 284.277 ns (0 allocations: 0 bytes)
# 155.669 ns (0 allocations: 0 bytes)
# 143.747 ns (0 allocations: 0 bytes)
# 462.041 μs (7 allocations: 20.42 KiB)
# 43.375 μs (0 allocations: 0 bytes)
# 37.541 μs (0 allocations: 0 bytes)
# 414.583 μs (0 allocations: 0 bytes)
# 414.500 μs (0 allocations: 0 bytes)
# 244.541 μs (0 allocations: 0 bytes)
# 210.250 μs (0 allocations: 0 bytes)
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

        GSEA.plot("", al, fe_, sc_, me_; la = Dict("title" => Dict("text" => strin(al))))

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
# 2.944 ms (13 allocations: 803.23 KiB)
# 2.676 ms (13 allocations: 803.23 KiB)
# 21.391 ms (13 allocations: 803.23 KiB)
# 21.517 ms (13 allocations: 803.23 KiB)
# 13.095 ms (13 allocations: 803.23 KiB)
# 11.457 ms (13 allocations: 803.23 KiB)
for al in AL_

    @test lastindex(GSEA.enrich(al, FE_, SO_, ME___)) === UE

    #@btime GSEA.enrich($al, FE_, SO_, ME___)

end

# ---- #

const SC = hcat(SO_, SO_ * 0.1, fill(0.8, lastindex(FE_)))

const RE = UE, size(SC, 2)

const OU = joinpath(homedir(), "Downloads")

# 9.548 ms (370 allocations: 5.51 MiB)
# 8.382 ms (370 allocations: 5.51 MiB)
# 52.180 ms (370 allocations: 5.51 MiB)
# 52.183 ms (370 allocations: 5.51 MiB)
# 30.802 ms (370 allocations: 5.51 MiB)
# 27.547 ms (370 allocations: 5.51 MiB)
#
# 9.355 ms (99 allocations: 6.03 MiB)
# 8.613 ms (99 allocations: 6.03 MiB)
# 64.771 ms (99 allocations: 6.03 MiB)
# 64.381 ms (99 allocations: 6.03 MiB)
# 39.669 ms (99 allocations: 6.03 MiB)
# 34.860 ms (99 allocations: 6.03 MiB)
for al in AL_

    en = GSEA.enrich(al, FE_, SC, ME___)

    @test size(en) === RE

    #@btime GSEA.enrich($al, FE_, SC, ME___)

    GSEA.write_plot(
        mkpath(joinpath(OU, strin(al))),
        FE_,
        SC,
        al,
        SE_,
        ME___,
        "Sample",
        ["Score", "Score x 0.1", "Constant"],
        en,
    )

end

# ---- #

const JS = joinpath(DA, "set_features.json")

# ---- #

# 2.148 ms (48968 allocations: 1.61 MiB)
# 16.805 ms (260901 allocations: 8.35 MiB)
# 2.172 ms (48974 allocations: 1.61 MiB)
for (fe_, mi, ma, fr, re) in (
    (String[], 33, 36, 0, 0),
    (unique!(vcat(values(Omics.Dic.rea(JS))...)), 33, 36, 0, 2),
    (["SHH", "XIST"], 1, 5656, 0, 2),
)

    se_, me___ = GSEA._read_set(JS, fe_, mi, ma, fr)

    @test lastindex(se_) === lastindex(me___) === re

    #@btime GSEA._read_set(JS, $fe_, $mi, $ma, $fr)

end

# ---- #

for (al, re) in zip(("ks", "ksa", "kliom", "kliop", "kli", "kli1"), AL_)

    @test GSEA._set_algorithm(al) == re

end

# ---- #

const TF = joinpath(DA, "feature_x_sample_x_number.tsv")

# ---- #

const OD = mkpath(joinpath(OU, "data_rank"))

GSEA.data_rank(OD, TF, JS)

const TD = Omics.Table.rea(joinpath(OD, "set_x_sample_x_enrichment.tsv"))

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

# 383.208 μs (2102 allocations: 2.39 MiB)
# 404.334 μs (2102 allocations: 2.39 MiB)
for al in (AL_[1], AL_[end])

    seed!(20231103)

    en_ = randn(100)

    ra = randn(100, 1000)

    GSEA._normalize_enrichment!(al, en_, ra)

    #@btime GSEA._normalize_enrichment!($al, $en_, $ra)

end

# ---- #

const FS = "set_x_statistic_x_result.tsv"

function test_result(ta, ue)

    @test size(ta, 1) === ue

end

function test_html(ou, uh)

    @test lastindex(filter!(re -> contains(re, r"html$"), readdir(ou))) === uh

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

const SU = Omics.Table.rea(joinpath(OS, FS))

test_result(SU, 50)

for (id, r1, r2, r3, r4) in (
    (1, "HALLMARK_PANCREAS_BETA_CELLS", -0.35266, -1.36616, 0.0200837),
    (2, "HALLMARK_PROTEIN_SECRETION", -0.272096, -1.25207, 0.0686192),
    (49, "HALLMARK_MYC_TARGETS_V1", 0.603356, 2.73998, 0.000262812),
    (50, "HALLMARK_MYC_TARGETS_V2", 0.866579, 3.36557, 0.000262812),
)

    @test SU[id, 1] === r1

    @test isapprox(SU[id, 2], r2; atol = 1e-6)

    @test isapprox(SU[id, 3], r3; atol = 1e-5)

    @test isapprox(SU[id, 4], r4; atol = 1e-7)

end

test_html(OS, 6)

# ---- #

const TT = joinpath(DA, "target_x_sample_x_number.tsv")

# ---- #

const OM = mkpath(joinpath(OU, "metric_rank"))

GSEA.metric_rank(OM, TT, TF, JS)

const TM = Omics.Table.rea(joinpath(OM, "feature_x_metric_x_score.tsv"))

@test size(TM) === (1000, 2)

@test names(TM) == ["Feature", "signal-to-noise-ratio"]

@test isapprox(TM[[1, 1000], 2], [1.8372355409610066, -1.7411005104346835])

test_result(Omics.Table.rea(joinpath(OM, FS)), 8)

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
