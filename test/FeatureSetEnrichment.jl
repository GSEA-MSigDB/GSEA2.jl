using Test: @test

using BioLab

# ---- #

const DA = joinpath(BioLab._DA, "FeatureSetEnrichment")

# ---- #

@test BioLab.Path.read(DA) ==
      ["c2.all.v7.1.symbols.gmt", "gene_x_statistic_x_number.tsv", "h.all.v7.1.symbols.gmt"]

# ---- #

const AL_ = (
    BioLab.FeatureSetEnrichment.KS(),
    BioLab.FeatureSetEnrichment.KSa(),
    BioLab.FeatureSetEnrichment.KLi(),
    BioLab.FeatureSetEnrichment.KLioM(),
    BioLab.FeatureSetEnrichment.KLioP(),
)

for (al, re) in zip(AL_, ("KS", "KSa", "KLi", "KLioM", "KLioP"))

    @test BioLab.FeatureSetEnrichment._make_string(al) == re

end

# ---- #

const SSC_ = [-2, -1, -0.5, 0, 0, 0.5, 1, 2, 3.4]

# ---- #

const ID = 1

for (ex, re) in (
    (-1, 0.5),
    (-1.0, 0.5),
    (0, 1),
    (0.0, 1),
    (1, 2),
    (1.0, 2),
    (2, 4),
    (2.0, 4),
    (3, 8),
    (3.0, 8),
    (4, 16),
    (4.0, 16),
    (0.1, 1.0717734625362931),
    (0.5, sqrt(2)),
)

    @test BioLab.FeatureSetEnrichment._index_absolute_exponentiate(SSC_, ID, ex) == re

    # 3.625 ns (0 allocations: 0 bytes)
    # 5.500 ns (0 allocations: 0 bytes)
    # 1.791 ns (0 allocations: 0 bytes)
    # 3.625 ns (0 allocations: 0 bytes)
    # 1.791 ns (0 allocations: 0 bytes)
    # 1.791 ns (0 allocations: 0 bytes)
    # 3.958 ns (0 allocations: 0 bytes)
    # 5.208 ns (0 allocations: 0 bytes)
    # 2.416 ns (0 allocations: 0 bytes)
    # 3.958 ns (0 allocations: 0 bytes)
    # 4.875 ns (0 allocations: 0 bytes)
    # 5.541 ns (0 allocations: 0 bytes)
    # 12.470 ns (0 allocations: 0 bytes)
    # 12.470 ns (0 allocations: 0 bytes)
    #@btime BioLab.FeatureSetEnrichment._index_absolute_exponentiate($SSC_, $ID, $ex)

end

# ---- #

const SIS_ = BitVector((true, false, true, false, true, true, false, false, true))

const N = length(SSC_)

# ---- #

for (ex, re) in (
    (1, (N, 4.0, 6.4)),
    (1.0, (N, 4.0, 6.4)),
    (2, (N, 4.0, 16.06)),
    (2.0, (N, 4.0, 16.06)),
    (0.1, (N, 4.0, 4.068020158853387)),
    (0.5, (N, 4.0, 4.672336016204768)),
)

    @test BioLab.FeatureSetEnrichment._sum_01(SSC_, ex, SIS_) == re

    # 9.843 ns (0 allocations: 0 bytes)
    # 10.010 ns (0 allocations: 0 bytes)
    # 21.314 ns (0 allocations: 0 bytes)
    # 34.659 ns (0 allocations: 0 bytes)
    # 57.249 ns (0 allocations: 0 bytes)
    # 57.206 ns (0 allocations: 0 bytes)
    #@btime BioLab.FeatureSetEnrichment._sum_01($SSC_, $ex, $SIS_)

end

# ---- #

for (ex, re) in (
    (1, (N, 10.4, 6.4)),
    (1.0, (N, 10.4, 6.4)),
    (2, (N, 22.06, 16.06)),
    (2.0, (N, 22.06, 16.06)),
    (0.1, (N, 7.13979362138968, 4.068020158853387)),
    (0.5, (N, 8.086549578577863, 4.672336016204768)),
)

    @test BioLab.FeatureSetEnrichment._sum_all1(SSC_, ex, SIS_) == re

    # 11.386 ns (0 allocations: 0 bytes)
    # 11.386 ns (0 allocations: 0 bytes)
    # 28.098 ns (0 allocations: 0 bytes)
    # 49.645 ns (0 allocations: 0 bytes)
    # 97.237 ns (0 allocations: 0 bytes)
    # 97.281 ns (0 allocations: 0 bytes)
    #@btime BioLab.FeatureSetEnrichment._sum_all1($SSC_, $ex, $SIS_)

end

# ---- #

const EX = 1

# ---- #

const CFE_ = ["K", "Q", "J", "X", "9", "8", "7", "6", "5", "4", "3", "2", "A"]

const CSC_ = [6.0, 5, 4, 3, 2, 1, 0, -1, -2, -3, -4, -5, -6]

const CFE1_ = ["A", "K"]

const CIS_ = in(Set(CFE1_)).(CFE_)

for (al, re) in zip(AL_, (-0.5, 0.0, 0.0, 0.0, 0.0))

    @test isapprox(
        BioLab.FeatureSetEnrichment._enrich!(al, CSC_, EX, CIS_, nothing),
        re;
        atol = 1e-15,
    )

    # 22.233 ns (0 allocations: 0 bytes)
    # 20.645 ns (0 allocations: 0 bytes)
    # 130.314 ns (0 allocations: 0 bytes)
    # 235.922 ns (0 allocations: 0 bytes)
    # 235.920 ns (0 allocations: 0 bytes)
    #@btime BioLab.FeatureSetEnrichment._enrich!($al, $CSC_, $EX, $CIS_, nothing)

end

# ---- #

for al in AL_

    BioLab.FeatureSetEnrichment.plot("", al, CFE_, CSC_, CFE1_, title_text = join((1:9..., 0))^10)

end

# ---- #

const MFE_, MSC_ =
    eachcol(BioLab.DataFrame.read(joinpath(DA, "gene_x_statistic_x_number.tsv"); select = [1, 2]))

reverse!(MFE_)

reverse!(MSC_)

const MFE1_ = BioLab.GMT.read(joinpath(DA, "c2.all.v7.1.symbols.gmt"))["COLLER_MYC_TARGETS_UP"]

const MIS_ = in(Set(MFE1_)).(MFE_)

const FE_X_SA_X_MSC = hcat(MSC_, MSC_ * 10.0, fill(0.8, length(MFE_)))

const MSE_FE1_ = BioLab.GMT.read(joinpath(DA, "h.all.v7.1.symbols.gmt"))

const MFE1___ = collect(values(MSE_FE1_))

# ---- #

for (al, re) in zip(
    AL_,
    (
        0.7651927829281453,
        0.41482514169516305,
        0.7736480596525319,
        0.7750661968892066,
        0.772229922415844,
    ),
)

    @test isapprox(
        BioLab.FeatureSetEnrichment._enrich!(al, MSC_, EX, MIS_, nothing),
        re;
        atol = 1e-12,
    )

    # 43.375 μs (0 allocations: 0 bytes)
    # 2.945 ms (108 allocations: 934.22 KiB)
    # 9.257 ms (358 allocations: 4.59 MiB)
    #
    # 37.166 μs (0 allocations: 0 bytes)
    # 2.640 ms (108 allocations: 934.22 KiB)
    # 8.343 ms (358 allocations: 4.59 MiB)
    #
    # 198.208 μs (0 allocations: 0 bytes)
    # 10.592 ms (108 allocations: 934.22 KiB)
    # 33.617 ms (358 allocations: 4.59 MiB)
    #
    # 349.667 μs (0 allocations: 0 bytes)
    # 18.327 ms (108 allocations: 934.22 KiB)
    # 55.971 ms (358 allocations: 4.59 MiB)
    #
    # 349.708 μs (0 allocations: 0 bytes)
    # 18.328 ms (108 allocations: 934.22 KiB)
    # 55.922 ms (358 allocations: 4.59 MiB)

    #@btime BioLab.FeatureSetEnrichment._enrich!($al, $MSC_, $EX, $MIS_, nothing)

    #@btime BioLab.FeatureSetEnrichment.enrich($al, $MFE_, $MSC_, $MFE1___)

    #@btime BioLab.FeatureSetEnrichment.enrich($al, $MFE_, $FE_X_SA_X_MSC, $MFE1___)

end

# ---- #

const MSE_ = collect(keys(MSE_FE1_))

const MSA_ = ["Score", "Score x 10", "Constant"]

const AL = BioLab.FeatureSetEnrichment.KS()

se_x_sa_x_en = BioLab.FeatureSetEnrichment.enrich(AL, MFE_, FE_X_SA_X_MSC, MFE1___)

@test BioLab.Error.@is_error BioLab.FeatureSetEnrichment.plot(
    "",
    AL,
    MFE_,
    FE_X_SA_X_MSC,
    MFE1___,
    "Sample",
    MSE_,
    MSA_,
    se_x_sa_x_en;
    ex = EX,
)

@test BioLab.FeatureSetEnrichment.plot(
    BioLab.TE,
    AL,
    MFE_,
    FE_X_SA_X_MSC,
    MFE1___,
    "Sample",
    MSE_,
    MSA_,
    se_x_sa_x_en;
    ex = EX,
) == BioLab.TE
