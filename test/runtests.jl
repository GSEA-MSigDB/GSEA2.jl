using Test: @test

using GSEA

# ----------------------------------------------------------------------------------------------- #

using Random: seed!

using Omics

# ---- #

const DD = pkgdir(GSEA, "data")

# ---- #

# 481.334 μs (7160 allocations: 473.16 KiB)
# 11.000 μs (116 allocations: 7.17 KiB)
# 10.833 μs (114 allocations: 6.84 KiB)

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

    va = Matrix(ta[:, 2:end])

    @test ta[:, 1][] === ph

    @test na_[1] === "Phenotype"

    @test all(startswith("Sample "), na_[2:end])

    @test eltype(va) === eltype(re)

    @test va[1, eachindex(re)] == re

    #@btime GSEA.read_cls($cl)

end

# ---- #

# 101.998 ms (71705 allocations: 23.67 MiB)

for (gc, re) in (("b.gct", (13321, 190)),)

    gc = joinpath(DD, gc)

    @test size(GSEA.read_gct(gc)) === re

    #@btime GSEA.read_gct($gc)

end

# ---- #

# 289.208 μs (7984 allocations: 1.12 MiB)
# 22.404 ms (537839 allocations: 62.61 MiB)

for (gm, re) in (("c.gmt", 50), ("d.gmt", 5529))

    gm = joinpath(DD, gm)

    se_me_ = GSEA.read_gmt(gm)

    @test typeof(se_me_) === Dict{String, Vector{String}}

    @test length(se_me_) === re

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

for (al, re) in zip(AL_, ("KS", "KSa", "KLioM", "KLioP", "KLi", "KLi1"))

    @test GSEA.strin(al) === re

    @test GSEA._set_algorithm(lowercase(re)) === al

end

# ---- #

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
# 69.629 ns (0 allocations: 0 bytes)
# 69.586 ns (0 allocations: 0 bytes)
# 8.299 ns (0 allocations: 0 bytes)
# 19.809 ns (0 allocations: 0 bytes)
# 115.695 ns (0 allocations: 0 bytes)
# 115.741 ns (0 allocations: 0 bytes)
# 8.291 ns (0 allocations: 0 bytes)
# 25.016 ns (0 allocations: 0 bytes)
# 78.832 ns (0 allocations: 0 bytes)
# 78.827 ns (0 allocations: 0 bytes)
# 10.761 ns (0 allocations: 0 bytes)
# 17.242 ns (0 allocations: 0 bytes)

const SG_ = [-2, -1, -0.5, 0, 0, 0.5, 1, 2, 3.4]

const IG_ = [true, false, true, false, true, true, false, false, true]

const N0 = -0.25

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

        @test GSEA._get_normalizer(al, SG_, ex, IG_) === re

        #@btime GSEA._get_normalizer($al, SG_, $ex, IG_)

    end

end

# ---- #

for (na, n1, re) in ((0.5, 1 / 3, -1.0),)

    @test GSEA._get_normalizer(na, n1) === re

end

# ---- #

const FM_, SM_ =
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
# 231.971 ns (6 allocations: 400 bytes)
# 17.285 ns (0 allocations: 0 bytes)
# 16.658 ns (0 allocations: 0 bytes)
# 284.155 ns (0 allocations: 0 bytes)
# 284.063 ns (0 allocations: 0 bytes)
# 155.717 ns (0 allocations: 0 bytes)
# 143.797 ns (0 allocations: 0 bytes)
# 460.334 μs (7 allocations: 20.42 KiB)
# 43.333 μs (0 allocations: 0 bytes)
# 37.500 μs (0 allocations: 0 bytes)
# 412.833 μs (0 allocations: 0 bytes)
# 414.041 μs (0 allocations: 0 bytes)
# 243.250 μs (0 allocations: 0 bytes)
# 210.291 μs (0 allocations: 0 bytes)

for (fe_, sc_, me_, re_) in (
    (
        ["K", "Q", "J", "X", "9", "8", "7", "6", "5", "4", "3", "2", "A"],
        [6.0, 5, 4, 3, 2, 1, 0, -1, -2, -3, -4, -5, -6],
        ["K", "A"],
        (-0.5, 0.0, 0.0, 0.0, 0.0, 0.0),
    ),
    (
        FM_,
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

@test GSEA._select_sort('a':'f', [1, NaN, 3, NaN, 5]) == (['e', 'c', 'a'], [5.0, 3.0, 1.0])

# ---- #

const SE_, ME___ = GSEA._separat(GSEA.read_gmt(joinpath(DD, "h.all.v7.1.symbols.gmt")))

# ---- #

# 3.022 ms (108 allocations: 934.22 KiB)
# 2.649 ms (108 allocations: 934.22 KiB)
# 17.202 ms (108 allocations: 934.22 KiB)
# 17.186 ms (108 allocations: 934.22 KiB)
# 10.127 ms (108 allocations: 934.22 KiB)
# 9.001 ms (108 allocations: 934.22 KiB)
#
# 3.051 ms (32 allocations: 1.86 MiB)
# 2.790 ms (32 allocations: 1.86 MiB)
# 21.514 ms (32 allocations: 1.86 MiB)
# 21.662 ms (32 allocations: 1.86 MiB)
# 13.287 ms (32 allocations: 1.86 MiB)
# 11.667 ms (32 allocations: 1.86 MiB)

for al in AL_

    en_ = GSEA.enrich(al, FM_, SM_, ME___)

    @test !issorted(en_)

    se_, en_ = GSEA._select_sort(SE_, en_)

    @test se_[1:2] == ["HALLMARK_MYC_TARGETS_V2", "HALLMARK_MYC_TARGETS_V1"]

    #@btime GSEA.enrich($al, FM_, SM_, ME___)

end

# ---- #

const DG = joinpath(tempdir(), "GSEA")

if isdir(DG)

    rm(DG; recursive = true)

end

mkdir(DG)

const FS = joinpath(DD, "set.json")

const FD = joinpath(DD, "data.tsv")

# ---- #

const OD = mkpath(joinpath(DG, "data_rank"))

const SD_, ED = GSEA.data_rank(OD, FD, FS; minimum_set_size = 15, maximum_set_size = 500)

@test isfile(joinpath(OD, "enrichment.tsv"))

@test SD_ == [
    "HALLMARK_ESTROGEN_RESPONSE_LATE",
    "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
    "HALLMARK_ESTROGEN_RESPONSE_EARLY",
    "HALLMARK_KRAS_SIGNALING_DN",
    "HALLMARK_IL2_STAT5_SIGNALING",
    "HALLMARK_APICAL_JUNCTION",
    "HALLMARK_HYPOXIA",
    "HALLMARK_GLYCOLYSIS",
]

@test findmax(ED) === (0.756249206577638, CartesianIndex(2, 2))

# ---- #

# 435.333 μs (1802 allocations: 3.03 MiB)
# 463.708 μs (1502 allocations: 3.69 MiB)

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

const OU = mkpath(joinpath(DG, "user_rank"))

GSEA.user_rank(
    OU,
    joinpath(DD, "metric.tsv"),
    FS;
    more_sets_to_plot = "HALLMARK_MYC_TARGETS_V1 HALLMARK_UV_RESPONSE_DN HALLMARK_UV_RESPONSE_UP ALIEN",
)

const RU = Omics.Table.rea(joinpath(OU, "result.tsv"))

test_result(RU, 50)

for (id, r1, r2, r3, r4) in (
    (45, "HALLMARK_PANCREAS_BETA_CELLS", -0.35266, -1.36616, 0.0200837),
    (33, "HALLMARK_PROTEIN_SECRETION", -0.272096, -1.25207, 0.0686192),
    (36, "HALLMARK_MYC_TARGETS_V1", 0.603356, 2.73998, 0.000262812),
    (10, "HALLMARK_MYC_TARGETS_V2", 0.866579, 3.36557, 0.000262812),
)

    @test RU[id, 1] === r1

    @test isapprox(RU[id, 2], r2; atol = 1e-6)

    # TODO: Investigate.

    #@test isapprox(RU[id, 3], r3; atol = 1e-5)

    #@test isapprox(RU[id, 4], r4; atol = 1e-7)

end

# ---- #

const OM = mkpath(joinpath(DG, "metric_rank"))

GSEA.metric_rank(
    OM,
    joinpath(DD, "target.tsv"),
    FD,
    FS;
    minimum_set_size = 15,
    maximum_set_size = 500,
)

const ME = Omics.Table.rea(joinpath(OM, "metric.tsv"))

@test size(ME) === (1000, 2)

@test names(ME) == ["Feature", "signal-to-noise-ratio"]

@test isapprox(sort(ME, 2)[[1, end], 2], [-1.8372355409610066, 1.7411005104346835])

const RU = Omics.Table.rea(joinpath(OM, "result.tsv"))

test_result(RU, 8)
