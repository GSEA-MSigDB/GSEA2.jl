using GSEA

using Test: @test

# ----------------------------------------------------------------------------------------------- #

using Omics

using CLSGCTGMT

# ---- #

const SC_ = [-2, -1, -0.5, 0, 0, 0.5, 1, 2, 3.4]

const IS_ = [true, false, true, false, true, true, false, false, true]

# ---- #

const N0 = -0.25

# 59.682 ns (0 allocations: 0 bytes)
# 59.681 ns (0 allocations: 0 bytes)
# 9.500 ns (0 allocations: 0 bytes)
# 21.314 ns (0 allocations: 0 bytes)
#
# 69.928 ns (0 allocations: 0 bytes)
# 69.970 ns (0 allocations: 0 bytes)
# 8.258 ns (0 allocations: 0 bytes)
# 19.831 ns (0 allocations: 0 bytes)
for (ex, re) in (
    (0.1, (N0, 0.24581982412836917)),
    (0.5, (N0, 0.21402570288861142)),
    (1, (N0, 0.15625)),
    (2, (N0, 0.06226650062266501)),
)

    @test GSEA._get_0_1(SC_, ex, IS_) === re

    @btime GSEA._get_0_1(SC_, $ex, IS_)

end

# ---- #

# 97.238 ns (0 allocations: 0 bytes)
# 97.281 ns (0 allocations: 0 bytes)
# 8.884 ns (0 allocations: 0 bytes)
# 28.098 ns (0 allocations: 0 bytes)
#
# 116.857 ns (0 allocations: 0 bytes)
# 116.858 ns (0 allocations: 0 bytes)
# 8.291 ns (0 allocations: 0 bytes)
# 25.016 ns (0 allocations: 0 bytes)
for (ex, re) in (
    (0.1, (0.14006007078470165, 0.24581982412836917)),
    (0.5, (0.12366213677204271, 0.21402570288861142)),
    (1, (0.09615384615384615, 0.15625)),
    (2, (0.04533091568449683, 0.06226650062266501)),
)

    @test GSEA._get_all_1(SC_, ex, IS_) === re

    @btime GSEA._get_all_1(SC_, $ex, IS_)

end

# ---- #

# 1.458 ns (0 allocations: 0 bytes)
#
# 2.083 ns (0 allocations: 0 bytes)
for (su, s1, re) in ((0.5, 1 / 3, -1.0),)

    @test GSEA._get_0(su, s1) === re

    @btime GSEA._get_0($su, $s1)

end

# ---- #

# 63.350 ns (0 allocations: 0 bytes)
# 63.350 ns (0 allocations: 0 bytes)
# 9.510 ns (0 allocations: 0 bytes)
# 20.687 ns (0 allocations: 0 bytes)
#
# 79.291 ns (0 allocations: 0 bytes)
# 79.296 ns (0 allocations: 0 bytes)
# 10.761 ns (0 allocations: 0 bytes)
# 17.285 ns (0 allocations: 0 bytes)
for (ex, re) in (
    (0.1, (0.24581982412836917)),
    (0.5, (0.21402570288861142)),
    (1, (0.15625)),
    (2, (0.06226650062266501)),
)

    @test GSEA._get_1(SC_, ex, IS_) === re

    @btime GSEA._get_1(SC_, $ex, IS_)

end

# ---- #

GSEA._clip

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
# 17.285 ns (0 allocations: 0 bytes)
# 16.658 ns (0 allocations: 0 bytes)
# 295.140 ns (0 allocations: 0 bytes)
# 295.295 ns (0 allocations: 0 bytes)
# 157.977 ns (0 allocations: 0 bytes)
# 147.877 ns (0 allocations: 0 bytes)
# 43.375 μs (0 allocations: 0 bytes)
# 37.541 μs (0 allocations: 0 bytes)
# 424.208 μs (0 allocations: 0 bytes)
# 424.167 μs (0 allocations: 0 bytes)
# 246.958 μs (0 allocations: 0 bytes)
# 215.250 μs (0 allocations: 0 bytes)
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

    is_ = map(in(Set(me_)), fe_)

    for (al, re) in zip(AL_, re_)

        ex = 1.0

        @test isapprox(GSEA._enrich!(al, sc_, ex, is_, nothing), re; atol = 1e-14)

        @btime GSEA._enrich!($al, $sc_, $ex, $is_, nothing)

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

GSEA._is_in!

# ---- #

const SE_ME_ = CLSGCTGMT.read_gmt(joinpath(DA, "h.all.v7.1.symbols.gmt"))

const SE_ = collect(keys(SE_ME_))

const ME___ = collect(values(SE_ME_))

const US = lastindex(SE_)

# ---- #

# 3.022 ms (108 allocations: 934.22 KiB)
# 2.649 ms (108 allocations: 934.22 KiB)
# 17.202 ms (108 allocations: 934.22 KiB)
# 17.186 ms (108 allocations: 934.22 KiB)
# 10.127 ms (108 allocations: 934.22 KiB)
# 9.001 ms (108 allocations: 934.22 KiB)
#
# 2.946 ms (13 allocations: 803.23 KiB)
# 2.683 ms (13 allocations: 803.23 KiB)
# 22.416 ms (13 allocations: 803.23 KiB)
# 22.511 ms (13 allocations: 803.23 KiB)
# 13.318 ms (13 allocations: 803.23 KiB)
# 11.693 ms (13 allocations: 803.23 KiB)
for al in AL_

    @test lastindex(GSEA.enrich(al, FE_, SO_, ME___)) === US

    @btime GSEA.enrich($al, FE_, SO_, ME___)

end

# ---- #

const SC = hcat(SO_, SO_ * 10.0, fill(0.8, lastindex(FE_)))

# ---- #

const RE = (US, size(SC, 2))

# 9.548 ms (370 allocations: 5.51 MiB)
# 8.382 ms (370 allocations: 5.51 MiB)
# 52.180 ms (370 allocations: 5.51 MiB)
# 52.183 ms (370 allocations: 5.51 MiB)
# 30.802 ms (370 allocations: 5.51 MiB)
# 27.547 ms (370 allocations: 5.51 MiB)
#
# 9.381 ms (99 allocations: 6.03 MiB)
# 8.559 ms (99 allocations: 6.03 MiB)
# 67.826 ms (99 allocations: 6.03 MiB)
# 67.611 ms (99 allocations: 6.03 MiB)
# 40.315 ms (99 allocations: 6.03 MiB)
# 35.746 ms (99 allocations: 6.03 MiB)
for al in AL_

    @test size(GSEA.enrich(al, FE_, SC, ME___)) === RE

    @btime GSEA.enrich($al, FE_, SC, ME___)

end
