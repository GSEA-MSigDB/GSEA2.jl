using GSEA

using Test: @test

# ----------------------------------------------------------------------------------------------- #

using Random: seed!

using Omics

using CLSGCTGMT

# ---- #

const FL_ = [-2, -1, -0.5, 0, 0, 0.5, 1, 2, 3.4]

# TODO: Benchmark Vector{Bool}
const BI_ = BitVector((1, 0, 1, 0, 1, 1, 0, 0, 1))

const RE = lastindex(FL_)

# ---- #

# 81.179 ns (0 allocations: 0 bytes)
# 81.190 ns (0 allocations: 0 bytes)
# 20.081 ns (0 allocations: 0 bytes)
# 19.998 ns (0 allocations: 0 bytes)
for (ex, re) in (
    (0.1, (RE, 0.24581982412836917)),
    (0.5, (RE, 0.21402570288861142)),
    (1, (RE, 0.15625)),
    (2, (RE, 0.06226650062266501)),
)

    @test GSEA._get_1_normalizer(FL_, ex, BI_) === re

    @btime GSEA._get_1_normalizer(FL_, $ex, BI_)

end

# ---- #

const RT = 0.25

# 71.954 ns (0 allocations: 0 bytes)
# 71.977 ns (0 allocations: 0 bytes)
# 19.747 ns (0 allocations: 0 bytes)
# 19.747 ns (0 allocations: 0 bytes)
for (ex, re) in (
    (0.1, (RE, RT, 0.24581982412836917)),
    (0.5, (RE, RT, 0.21402570288861142)),
    (1, (RE, RT, 0.15625)),
    (2, (RE, RT, 0.06226650062266501)),
)

    @test GSEA._get_0_1_normalizer(FL_, ex, BI_) === re

    @btime GSEA._get_0_1_normalizer(FL_, $ex, BI_)

end

# ---- #

# 119.054 ns (0 allocations: 0 bytes)
# 119.009 ns (0 allocations: 0 bytes)
# 25.602 ns (0 allocations: 0 bytes)
# 26.857 ns (0 allocations: 0 bytes)
for (ex, re) in (
    (0.1, (RE, 0.14006007078470165, 0.24581982412836917)),
    (0.5, (RE, 0.12366213677204271, 0.21402570288861142)),
    (1, (RE, 0.09615384615384615, 0.15625)),
    (2, (RE, 0.04533091568449683, 0.06226650062266501)),
)

    @test GSEA._get_all_1_normalizer(FL_, ex, BI_) === re

    @btime GSEA._get_all_1_normalizer(FL_, $ex, BI_)

end

# ---- #

# 2.084 ns (0 allocations: 0 bytes)
for (su, s1, re) in ((0.5, 1 / 3, -1.0),)

    @test GSEA._get_0_normalizer(su, s1) === re

    @btime GSEA._get_0_normalizer($su, $s1)

end

# ---- #

const AL_ = GSEA.KS(), GSEA.KSa(), GSEA.KLi1(), GSEA.KLi(), GSEA.KLioM(), GSEA.KLioP()

# ---- #

const DA = pkgdir(GSEA, "data")

const FE_, SC_ = eachcol(
    reverse!(
        Omics.Table.rea(joinpath(DA, "gene_x_statistic_x_number.tsv"); select = [1, 2]),
    ),
)

# ---- #

# 31.942 ns (0 allocations: 0 bytes)
# 31.816 ns (0 allocations: 0 bytes)
# 165.308 ns (0 allocations: 0 bytes)
# 258.069 ns (0 allocations: 0 bytes)
# 403.750 ns (0 allocations: 0 bytes)
# 403.750 ns (0 allocations: 0 bytes)
# 45.208 μs (0 allocations: 0 bytes)
# 37.500 μs (0 allocations: 0 bytes)
# 164.833 μs (0 allocations: 0 bytes)
# 186.208 μs (0 allocations: 0 bytes)
# 325.833 μs (0 allocations: 0 bytes)
# 326.042 μs (0 allocations: 0 bytes)
for (fe_, sc_, me_, re_) in (
    (
        ["K", "Q", "J", "X", "9", "8", "7", "6", "5", "4", "3", "2", "A"],
        [6.0, 5, 4, 3, 2, 1, 0, -1, -2, -3, -4, -5, -6],
        ["K", "A"],
        (-0.5, 0.0, 0.0, 0.0, 0.0, 0.0),
    ),
    (
        FE_,
        SC_,
        CLSGCTGMT.read_gmt(joinpath(DA, "c2.all.v7.1.symbols.gmt"))["COLLER_MYC_TARGETS_UP"],
        (
            0.7651927829281453,
            0.41482514169516305,
            #0.8524266036047564,
            #0.7736480596525319,
            #0.7750661968892066,
            #0.772229922415844,
        ),
    ),
)

    is_ = map(in(Set(me_)), fe_)

    for (al, re) in zip(AL_, re_)

        ex = 1.0

        @test isapprox(GSEA._enrich!(al, sc_, ex, is_, nothing), re; atol = 1e-15)

        GSEA.plot(
            "",
            al,
            fe_,
            sc_,
            me_;
            la = Dict("title" => Dict("text" => string(al)[6:(end - 2)])),
        )

        @btime GSEA._enrich!($al, $sc_, $ex, $is_, nothing)

    end

end

# ---- #

const SE_ME_ = CLSGCTGMT.read_gmt(joinpath(DA, "h.all.v7.1.symbols.gmt"))

const SE_ = collect(keys(SE_ME_))

const ME___ = collect(values(SE_ME_))

# ---- #

# 3.022 ms (108 allocations: 934.22 KiB)
# 2.649 ms (108 allocations: 934.22 KiB)
# 9.001 ms (108 allocations: 934.22 KiB)
# 10.127 ms (108 allocations: 934.22 KiB)
# 17.202 ms (108 allocations: 934.22 KiB)
# 17.186 ms (108 allocations: 934.22 KiB)
for al in AL_

    @btime GSEA.enrich($al, FE_, SC_, ME___)

end

# ---- #

const SC = hcat(SC_, SC_ * 10.0, fill(0.8, lastindex(FE_)))

# ---- #

# 9.548 ms (370 allocations: 5.51 MiB)
# 8.382 ms (370 allocations: 5.51 MiB)
# 27.547 ms (370 allocations: 5.51 MiB)
# 30.802 ms (370 allocations: 5.51 MiB)
# 52.180 ms (370 allocations: 5.51 MiB)
# 52.183 ms (370 allocations: 5.51 MiB)
for al in AL_

    @btime GSEA.enrich($al, FE_, SC, ME___)

end

# ---- #

const AL = GSEA.KS()

GSEA.plot(
    tempdir(),
    AL,
    FE_,
    SC,
    ME___,
    "Sample",
    SE_,
    ["Score", "Score x 10", "Constant"],
    GSEA.enrich(AL, FE_, SC, ME___),
)
