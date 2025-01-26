using Test: @test

using GSEA

# ----------------------------------------------------------------------------------------------- #

using Random: seed!

using Omics

# ---- #

const DD = pkgdir(GSEA, "data")

# ---- #

# 401.583 μs (7160 allocations: 473.16 KiB)
# 10.458 μs (116 allocations: 7.17 KiB)
# 10.166 μs (114 allocations: 6.84 KiB)
for (cl, ph, re) in (
    (
        "CCLE_mRNA_20Q2_no_haem_phen.cls",
        "HER2",
        [1.087973, -1.358492, -1.178614, -0.77898, 0.157222, 1.168224, -0.360195, 0.608629],
    ),
    ("GSE76137.cls", "Proliferating_Arrested", [1, 2, 1, 2, 1, 2]),
    ("a.cls", "CNTRL_LPS", [1, 1, 1, 2, 2, 2]),
)

    cl = joinpath(DD, cl)

    ta = GSEA.read_cls(cl)

    na_ = names(ta)

    ma = Matrix(ta[:, 2:end])

    @test ta[:, 1][] === ph

    @test na_[1] === "Phenotype"

    @test all(startswith("Sample "), na_[2:end])

    @test eltype(ma) === eltype(re)

    @test ma[1, eachindex(re)] == re

    #@btime GSEA.read_cls($cl)

end

# ---- #

# 102.867 ms (71705 allocations: 23.67 MiB)
for (gc, re) in (("b.gct", (13321, 190)),)

    gc = joinpath(DD, gc)

    @test size(GSEA.read_gct(gc)) === re

    #@btime GSEA.read_gct($gc)

end

# ---- #

# 286.000 μs (7984 allocations: 1.12 MiB)
# 22.167 ms (537839 allocations: 62.61 MiB)
for (gm, re) in (("c.gmt", 50), ("d.gmt", 5529))

    gm = joinpath(DD, gm)

    se_ge_ = GSEA.read_gmt(gm)

    @test typeof(se_ge_) === Dict{String, Vector{String}}

    @test length(se_ge_) === re

    #@btime GSEA.read_gmt($gm)

end

# ---- #

GSEA.cls

# ---- #

GSEA.gct

# ---- #

GSEA.gmt

# ---- #

const AL_ = GSEA.KS(), GSEA.KSa(), GSEA.KLioM(), GSEA.KLioP(), GSEA.KLi(), GSEA.KLi1()

# ---- #

for (al, re) in zip(AL_, ("KS", "KSa", "KLioM", "KLioP", "KLi", "KLi1"))

    @test GSEA.strin(al) == re

    @test GSEA._set_algorithm(lowercase(re)) == al

end

# ---- #

const SG_ = [-2, -1, -0.5, 0, 0, 0.5, 1, 2, 3.4]

const II_ = [true, false, true, false, true, true, false, false, true]

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
# 8.250 ns (0 allocations: 0 bytes)
# 19.851 ns (0 allocations: 0 bytes)
# 116.884 ns (0 allocations: 0 bytes)
# 116.857 ns (0 allocations: 0 bytes)
# 8.299 ns (0 allocations: 0 bytes)
# 25.016 ns (0 allocations: 0 bytes)
# 79.260 ns (0 allocations: 0 bytes)
# 79.291 ns (0 allocations: 0 bytes)
# 10.750 ns (0 allocations: 0 bytes)
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

        @test GSEA._get_normalizer(al, SG_, ex, II_) === re

        #@btime GSEA._get_delta($al, SG_, $ex, II_)

    end

end

# ---- #

# 1.458 ns (0 allocations: 0 bytes)
#
# 2.083 ns (0 allocations: 0 bytes)
for (na, n1, re) in ((0.5, 1 / 3, -1.0),)

    @test GSEA._get_normalizer(na, n1) === re

    #@btime GSEA._get_delta($na, $n1)

end

# ---- #

const FE_, SM_ =
    eachcol(reverse!(Omics.Table.rea(joinpath(DD, "myc.tsv"); select = [1, 2])))

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
# 233.424 ns (6 allocations: 400 bytes)
# 17.285 ns (0 allocations: 0 bytes)
# 16.658 ns (0 allocations: 0 bytes)
# 284.452 ns (0 allocations: 0 bytes)
# 284.574 ns (0 allocations: 0 bytes)
# 155.978 ns (0 allocations: 0 bytes)
# 143.796 ns (0 allocations: 0 bytes)
# 464.958 μs (7 allocations: 20.42 KiB)
# 43.375 μs (0 allocations: 0 bytes)
# 37.500 μs (0 allocations: 0 bytes)
# 413.041 μs (0 allocations: 0 bytes)
# 413.000 μs (0 allocations: 0 bytes)
# 243.209 μs (0 allocations: 0 bytes)
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
        SM_,
        GSEA.read_gmt(joinpath(DD, "c2.all.v7.1.symbols.gmt"))["COLLER_MYC_TARGETS_UP"],
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

    ii_ = GSEA._is_in(fe_, me_)

    @test typeof(ii_) === Vector{Bool}

    #@btime GSEA._is_in($fe_, $me_)

    for (al, re) in zip(AL_, re_)

        ex = 1.0

        @test isapprox(GSEA._enrich!(al, sc_, ex, ii_, nothing), re; atol = 1e-11)

        #@btime GSEA._enrich!($al, $sc_, $ex, $ii_, nothing)

        GSEA.plot(
            "",
            al,
            fe_,
            sc_,
            me_;
            la = Dict("title" => Dict("text" => GSEA.strin(al))),
        )

    end

end

# ---- #

const SE_ME_ = GSEA.read_gmt(joinpath(DD, "h.all.v7.1.symbols.gmt"))

const SE_, ME___ = GSEA._separat(SE_ME_)

const US = lastindex(SE_)

# ---- #

@test GSEA._select_sort('a':'f', [1, NaN, 3, NaN, 5]) == (['e', 'c', 'a'], [5.0, 3.0, 1.0])

# ---- #

# 3.022 ms (108 allocations: 934.22 KiB)
# 2.649 ms (108 allocations: 934.22 KiB)
# 17.202 ms (108 allocations: 934.22 KiB)
# 17.186 ms (108 allocations: 934.22 KiB)
# 10.127 ms (108 allocations: 934.22 KiB)
# 9.001 ms (108 allocations: 934.22 KiB)
#
# 2.939 ms (13 allocations: 803.23 KiB)
# 2.681 ms (13 allocations: 803.23 KiB)
# 21.308 ms (13 allocations: 803.23 KiB)
# 21.280 ms (13 allocations: 803.23 KiB)
# 13.064 ms (13 allocations: 803.23 KiB)
# 11.444 ms (13 allocations: 803.23 KiB)
for al in AL_

    @test lastindex(GSEA.enrich(al, FE_, SM_, ME___)) === US

    #@btime GSEA.enrich($al, FE_, SM_, ME___)

end

# ---- #

const SC = hcat(SM_, SM_ * 0.1, fill(0.8, lastindex(FE_)))

const RE = US, size(SC, 2)

const DU = joinpath(homedir(), "Downloads", "GSEA")

# 9.548 ms (370 allocations: 5.51 MiB)
# 8.382 ms (370 allocations: 5.51 MiB)
# 52.180 ms (370 allocations: 5.51 MiB)
# 52.183 ms (370 allocations: 5.51 MiB)
# 30.802 ms (370 allocations: 5.51 MiB)
# 27.547 ms (370 allocations: 5.51 MiB)
#
# 9.386 ms (99 allocations: 6.03 MiB)
# 8.594 ms (99 allocations: 6.03 MiB)
# 64.708 ms (99 allocations: 6.03 MiB)
# 64.689 ms (99 allocations: 6.03 MiB)
# 39.662 ms (99 allocations: 6.03 MiB)
# 34.834 ms (99 allocations: 6.03 MiB)
for al in AL_

    en = stack(map(sc_ -> GSEA.enrich(al, FE_, sc_, ME___; fr = 0.0), eachcol(SC)))

    @test size(en) === RE

    #@btime GSEA.enrich($al, FE_, SC, ME___)

    GSEA.write_plot(
        mkpath(joinpath(DU, GSEA.strin(al))),
        al,
        FE_,
        SC,
        SE_,
        ME___,
        "Sample",
        ["Score", "Score x 0.1", "Constant"],
        en;
        up = 8,
    )

end

# ---- #

const JS = joinpath(DD, "set.json")

const TD = joinpath(DD, "data.tsv")

# ---- #

const OD = mkpath(joinpath(DU, "data_rank"))

GSEA.data_rank(OD, TD, JS)

const TE = Omics.Table.rea(joinpath(OD, "enrichment.tsv"))

@test size(TE) === (50, 10)

#@test TE[!, "Set"] == [
#    "HALLMARK_ESTROGEN_RESPONSE_LATE",
#    "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
#    "HALLMARK_ESTROGEN_RESPONSE_EARLY",
#    "HALLMARK_KRAS_SIGNALING_DN",
#    "HALLMARK_IL2_STAT5_SIGNALING",
#    "HALLMARK_APICAL_JUNCTION",
#    "HALLMARK_HYPOXIA",
#    "HALLMARK_GLYCOLYSIS",
#]

# ---- #

# 490.541 μs (2402 allocations: 3.17 MiB)
# 399.750 μs (2002 allocations: 3.15 MiB)
for al in (AL_[1], AL_[end])

    seed!(20231103)

    en_ = randn(100)

    ra = randn(100, 1000)

    GSEA._normalize_enrichment!(al, en_, ra)

    #@btime GSEA._normalize_enrichment!($al, $en_, $ra)

end

# ---- #

function test_result(ta, US)

    @test size(ta, 1) === US

end

# ---- #

const OS = mkpath(joinpath(DU, "user_rank"))

GSEA.user_rank(
    OS,
    joinpath(DD, "metric.tsv"),
    JS;
    number_of_sets_to_plot = 8,
    more_sets_to_plot = "HALLMARK_MYC_TARGETS_V1 HALLMARK_UV_RESPONSE_DN HALLMARK_UV_RESPONSE_UP ALIEN",
)

const SU = Omics.Table.rea(joinpath(OS, "result.tsv"))

test_result(SU, 50)

for (id, r1, r2, r3, r4) in (
    (1, "HALLMARK_PANCREAS_BETA_CELLS", -0.35266, -1.36616, 0.0200837),
    (2, "HALLMARK_PROTEIN_SECRETION", -0.272096, -1.25207, 0.0686192),
    (49, "HALLMARK_MYC_TARGETS_V1", 0.603356, 2.73998, 0.000262812),
    (50, "HALLMARK_MYC_TARGETS_V2", 0.866579, 3.36557, 0.000262812),
)

    @test SU[id, 1] === r1

    @test isapprox(SU[id, 2], r2; atol = 1e-6)

    #@test isapprox(SU[id, 3], r3; atol = 1e-5)

    #@test isapprox(SU[id, 4], r4; atol = 1e-7)

end

# ---- #

const TT = joinpath(DD, "target.tsv")

const OM = mkpath(joinpath(DU, "metric_rank"))

GSEA.metric_rank(OM, TT, TD, JS; number_of_sets_to_plot = 8)

const TM = Omics.Table.rea(joinpath(OM, "metric.tsv"))

@test size(TM) === (1000, 2)

@test names(TM) == ["Feature", "signal-to-noise-ratio"]

@test isapprox(sort(TM, 2)[[1, 1000], 2], [-1.8372355409610066, 1.7411005104346835])

test_result(Omics.Table.rea(joinpath(OM, "result.tsv")), 50)
