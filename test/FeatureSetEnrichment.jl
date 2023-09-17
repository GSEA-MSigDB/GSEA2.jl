using Test: @test

using BioLab

using GSEA: FeatureSetEnrichment

# ---- #

const DA = joinpath(dirname(@__DIR__), "data", "FeatureSetEnrichment")

@test BioLab.Path.read(DA) == [
    "c2.all.v7.1.symbols.gmt",
    "coller_auc.tsv",
    "gene_x_statistic_x_number.tsv",
    "h.all.v2022.1.Hs.symbols.gmt",
    "h.all.v7.1.symbols.gmt",
]

# ---- #

const AL_ = (
    FeatureSetEnrichment.KS(),
    FeatureSetEnrichment.KSa(),
    FeatureSetEnrichment.KLi1(),
    FeatureSetEnrichment.KLi(),
    FeatureSetEnrichment.KLioM(),
    FeatureSetEnrichment.KLioP(),
)

# ---- #

for (al, re) in zip(AL_, ("KS", "KSa", "KLi1", "KLi", "KLioM", "KLioP"))

    @test FeatureSetEnrichment._make_string(al) == re

end

# ---- #

const SSC_ = [-2, -1, -0.5, 0, 0, 0.5, 1, 2, 3.4]

# ---- #

const ID = 1

for (ex, re) in (
    (-1, 0.5),
    (-1.0, 0.5),
    (0, 1.0),
    (0.0, 1.0),
    (1, 2.0),
    (1.0, 2.0),
    (2, 4.0),
    (2.0, 4.0),
    (3, 8.0),
    (3.0, 8.0),
    (4, 16.0),
    (4.0, 16.0),
    (0.1, 1.0717734625362931),
    (0.5, sqrt(2)),
)

    @test FeatureSetEnrichment._index_absolute_exponentiate(SSC_, ID, ex) === re

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
    #@btime FeatureSetEnrichment._index_absolute_exponentiate($SSC_, $ID, $ex)

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

    @test FeatureSetEnrichment._sum_01(SSC_, ex, SIS_) === re

    # 9.843 ns (0 allocations: 0 bytes)
    # 7.375 ns (0 allocations: 0 bytes)
    # 21.063 ns (0 allocations: 0 bytes)
    # 28.811 ns (0 allocations: 0 bytes)
    # 56.021 ns (0 allocations: 0 bytes)
    # 56.063 ns (0 allocations: 0 bytes)
    #@btime FeatureSetEnrichment._sum_01($SSC_, $ex, $SIS_)

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

    @test FeatureSetEnrichment._sum_all1(SSC_, ex, SIS_) === re

    # 11.052 ns (0 allocations: 0 bytes)
    # 10.511 ns (0 allocations: 0 bytes)
    # 28.098 ns (0 allocations: 0 bytes)
    # 49.047 ns (0 allocations: 0 bytes)
    # 96.636 ns (0 allocations: 0 bytes)
    # 96.724 ns (0 allocations: 0 bytes)
    #@btime FeatureSetEnrichment._sum_all1($SSC_, $ex, $SIS_)

end

# ---- #

const EX = 1

# ---- #

const CFE_ = ["K", "Q", "J", "X", "9", "8", "7", "6", "5", "4", "3", "2", "A"]

const CSC_ = [6.0, 5, 4, 3, 2, 1, 0, -1, -2, -3, -4, -5, -6]

const CFE1_ = ["K", "A"]

const CIS_ = in(Set(CFE1_)).(CFE_)

# ---- #

for (al, re) in zip(AL_, (-0.5, 0, 0, 0, 0, 0))

    @test isapprox(FeatureSetEnrichment._enrich!(al, CSC_, EX, CIS_, nothing), re; atol = 1e-15)

    # 19.433 ns (0 allocations: 0 bytes)
    # 17.577 ns (0 allocations: 0 bytes)
    # 153.895 ns (0 allocations: 0 bytes)
    # 129.467 ns (0 allocations: 0 bytes)
    # 241.137 ns (0 allocations: 0 bytes)
    # 235.118 ns (0 allocations: 0 bytes)
    #@btime FeatureSetEnrichment._enrich!($al, $CSC_, $EX, $CIS_, nothing)

end

# ---- #

for al in AL_

    FeatureSetEnrichment.plot("", al, CFE_, CSC_, CFE1_; title_text = join((1:9..., 0))^10)

end

# ---- #

const PFE_, PSC_ = eachcol(BioLab.DataFrame.read(joinpath(DA, "Coller_auc.tsv"); select = [1, 2]))

const PIS_ =
    in(
        Set(
            BioLab.GMT.read(joinpath(DA, "h.all.v2022.1.Hs.symbols.gmt"))["HALLMARK_MYC_TARGETS_V1"],
        ),
    ).(PFE_)

for (al, re) in zip(AL_, (
    0.6823,
    0.3988,
    0.817,
    0.7171,
    0.7287,
    0.7055,
    #0.3643,
    #0.3528,
))

    @test round(FeatureSetEnrichment._enrich!(al, PSC_, EX, PIS_, nothing); digits = 4) === re

    # 43.250 μs (0 allocations: 0 bytes)
    # 37.083 μs (0 allocations: 0 bytes)
    # 193.667 μs (0 allocations: 0 bytes)
    # 196.708 μs (0 allocations: 0 bytes)
    # 357.459 μs (0 allocations: 0 bytes)
    # 348.291 μs (0 allocations: 0 bytes)
    #@btime FeatureSetEnrichment._enrich!($al, $PSC_, $EX, $PIS_, nothing)

end

# ---- #

const MFE_, MSC_ =
    eachcol(BioLab.DataFrame.read(joinpath(DA, "gene_x_statistic_x_number.tsv"); select = [1, 2]))

reverse!(MFE_)

reverse!(MSC_)

const MIS_ =
    in(
        Set(BioLab.GMT.read(joinpath(DA, "c2.all.v7.1.symbols.gmt"))["COLLER_MYC_TARGETS_UP"]),
    ).(MFE_)

const FE_X_SA_X_MSC = hcat(MSC_, MSC_ * 10.0, fill(0.8, length(MFE_)))

const MSE_FE1_ = BioLab.GMT.read(joinpath(DA, "h.all.v7.1.symbols.gmt"))

const MFE1___ = collect(values(MSE_FE1_))

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
        #0.3875330984446027,
        #0.38611496120792593,
    ),
)

    @test isapprox(FeatureSetEnrichment._enrich!(al, MSC_, EX, MIS_, nothing), re; atol = 1e-12)

    # 43.375 μs (0 allocations: 0 bytes)
    # 37.166 μs (0 allocations: 0 bytes)
    # 193.584 μs (0 allocations: 0 bytes)
    # 196.375 μs (0 allocations: 0 bytes)
    # 348.750 μs (0 allocations: 0 bytes)
    # 348.458 μs (0 allocations: 0 bytes)
    #@btime FeatureSetEnrichment._enrich!($al, $MSC_, $EX, $MIS_, nothing)

end

# ---- #

for al in AL_

    # 2.943 ms (108 allocations: 934.22 KiB)
    # 2.633 ms (108 allocations: 934.22 KiB)
    # 10.407 ms (108 allocations: 934.22 KiB)
    # 10.477 ms (108 allocations: 934.22 KiB)
    # 18.656 ms (108 allocations: 934.22 KiB)
    # 18.227 ms (108 allocations: 934.22 KiB)
    #@btime FeatureSetEnrichment.enrich($al, $MFE_, $MSC_, $MFE1___)

end

# ---- #

for al in AL_

    # 9.231 ms (358 allocations: 4.59 MiB)
    # 8.330 ms (358 allocations: 4.59 MiB)
    # 32.917 ms (358 allocations: 4.59 MiB)
    # 33.163 ms (358 allocations: 4.59 MiB)
    # 55.555 ms (358 allocations: 4.59 MiB)
    # 55.553 ms (358 allocations: 4.59 MiB)
    #@btime FeatureSetEnrichment.enrich($al, $MFE_, $FE_X_SA_X_MSC, $MFE1___)

end

# ---- #

const MSE_ = collect(keys(MSE_FE1_))

const MSA_ = ["Score", "Score x 10", "Constant"]

const AL = FeatureSetEnrichment.KS()

const SE_X_SA_X_EN = FeatureSetEnrichment.enrich(AL, MFE_, FE_X_SA_X_MSC, MFE1___)

# ---- #

@test BioLab.Error.@is FeatureSetEnrichment.plot(
    "",
    AL,
    MFE_,
    FE_X_SA_X_MSC,
    MFE1___,
    "Sample",
    MSE_,
    MSA_,
    SE_X_SA_X_EN;
    ex = EX,
)

# ---- #

@test FeatureSetEnrichment.plot(
    BioLab.TE,
    AL,
    MFE_,
    FE_X_SA_X_MSC,
    MFE1___,
    "Sample",
    MSE_,
    MSA_,
    SE_X_SA_X_EN;
    ex = EX,
) === BioLab.TE
