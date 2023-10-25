using Test: @test

using Nucleus

using GSEA: FeatureSetEnrichment

# ---- #

const DA = joinpath(dirname(@__DIR__), "data", "FeatureSetEnrichment")

# ---- #

@test Nucleus.Path.read(DA) == [
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

    # 549.021 ns (11 allocations: 424 bytes)
    # 559.360 ns (11 allocations: 424 bytes)
    # 568.158 ns (11 allocations: 424 bytes)
    # 555.108 ns (11 allocations: 424 bytes)
    # 566.799 ns (11 allocations: 424 bytes)
    # 560.584 ns (11 allocations: 424 bytes)
    #@btime FeatureSetEnrichment._make_string($al)

end

# ---- #

const SSC_ = [-2, -1, -0.5, 0, 0, 0.5, 1, 2, 3.4]

# ---- #

const SC = SSC_[1]

# ---- #

for ex in (-1, 0, 1, 2, 3, 4, 0.1, 0.5)

    @test FeatureSetEnrichment._absolute_exponentiate(SC, ex) === abs(SC)^ex

    # 3.666 ns (0 allocations: 0 bytes)
    # 1.792 ns (0 allocations: 0 bytes)
    # 1.500 ns (0 allocations: 0 bytes)
    # 3.959 ns (0 allocations: 0 bytes)
    # 2.375 ns (0 allocations: 0 bytes)
    # 4.875 ns (0 allocations: 0 bytes)
    # 11.386 ns (0 allocations: 0 bytes)
    # 11.386 ns (0 allocations: 0 bytes)
    #@btime FeatureSetEnrichment._absolute_exponentiate(SC, $ex)

end

# ---- #

const SIS_ = BitVector((true, false, true, false, true, true, false, false, true))

# ---- #

const N = lastindex(SSC_)

# ---- #

for (id, fu) in enumerate((
    FeatureSetEnrichment._get_0_1_normalizer,
    FeatureSetEnrichment._get_1_normalizer,
    FeatureSetEnrichment._get_all_1_normalizer,
))

    for (ex, re_...) in (
        (1, (0.25, 0.15625), (0.15625,), (0.09615384615384615, 0.15625)),
        (
            2,
            (0.25, 0.06226650062266501),
            (0.06226650062266501,),
            (0.04533091568449683, 0.06226650062266501),
        ),
        (
            0.1,
            (0.25, 0.24581982412836917),
            (0.24581982412836917,),
            (0.14006007078470165, 0.24581982412836917),
        ),
        (
            0.5,
            (0.25, 0.21402570288861142),
            (0.21402570288861142,),
            (0.12366213677204271, 0.21402570288861142),
        ),
    )

        @test fu(SSC_, ex, SIS_) === (N, re_[id]...)

        # 9.510 ns (0 allocations: 0 bytes)
        # 21.314 ns (0 allocations: 0 bytes)
        # 59.700 ns (0 allocations: 0 bytes)
        # 59.700 ns (0 allocations: 0 bytes)
        # 9.510 ns (0 allocations: 0 bytes)
        # 20.687 ns (0 allocations: 0 bytes)
        # 63.350 ns (0 allocations: 0 bytes)
        # 63.350 ns (0 allocations: 0 bytes)
        # 8.884 ns (0 allocations: 0 bytes)
        # 28.098 ns (0 allocations: 0 bytes)
        # 97.281 ns (0 allocations: 0 bytes)
        # 97.324 ns (0 allocations: 0 bytes)
        #@btime $fu(SSC_, $ex, SIS_)

    end

end

# ---- #

@test FeatureSetEnrichment._get_0_normalizer(1 / 2, 1 / 3) === -1.0

# ---- #

# 0.875 ns (0 allocations: 0 bytes)
#@btime FeatureSetEnrichment._get_0_normalizer(1 / 2, 1 / 3)

# ---- #

const EX = 1

# ---- #

const CFE_ = ["K", "Q", "J", "X", "9", "8", "7", "6", "5", "4", "3", "2", "A"]

# ---- #

const CSC_ = [6.0, 5, 4, 3, 2, 1, 0, -1, -2, -3, -4, -5, -6]

# ---- #

const CFE1_ = ["K", "A"]

# ---- #

const CIS_ = in(Set(CFE1_)).(CFE_)

# ---- #

for (al, re) in zip(AL_, (-0.5, 0, 0, 0, 0, 0))

    @test isapprox(FeatureSetEnrichment._enrich!(al, CSC_, EX, CIS_, nothing), re; atol = 1e-15)

    # 19.455 ns (0 allocations: 0 bytes)
    # 17.869 ns (0 allocations: 0 bytes)
    # 117.569 ns (0 allocations: 0 bytes)
    # 126.717 ns (0 allocations: 0 bytes)
    # 225.789 ns (0 allocations: 0 bytes)
    # 225.873 ns (0 allocations: 0 bytes)
    #@btime FeatureSetEnrichment._enrich!($al, CSC_, EX, CIS_, nothing)

end

# ---- #

for al in AL_

    FeatureSetEnrichment.plot("", al, CFE_, CSC_, CFE1_; title_text = join((1:9..., 0))^10)

end

# ---- #

const PFE_, PSC_ = eachcol(Nucleus.DataFrame.read(joinpath(DA, "Coller_auc.tsv"); select = [1, 2]))

# ---- #

const PIS_ =
    in(
        Set(
            Nucleus.GMT.read(joinpath(DA, "h.all.v2022.1.Hs.symbols.gmt"))["HALLMARK_MYC_TARGETS_V1"],
        ),
    ).(PFE_)

# ---- #

for (al, re) in zip(AL_, (0.6823, 0.3988, 0.817, 0.7171, 0.7287, 0.7055))

    @test round(FeatureSetEnrichment._enrich!(al, PSC_, EX, PIS_, nothing); digits = 4) === re

    # 43.208 μs (0 allocations: 0 bytes)
    # 37.541 μs (0 allocations: 0 bytes)
    # 165.125 μs (0 allocations: 0 bytes)
    # 186.584 μs (0 allocations: 0 bytes)
    # 326.625 μs (0 allocations: 0 bytes)
    # 326.042 μs (0 allocations: 0 bytes)
    #@btime FeatureSetEnrichment._enrich!($al, PSC_, EX, PIS_, nothing)

end

# ---- #

const MFE_, MSC_ =
    eachcol(Nucleus.DataFrame.read(joinpath(DA, "gene_x_statistic_x_number.tsv"); select = [1, 2]))

# ---- #

reverse!(MFE_)

# ---- #

reverse!(MSC_)

# ---- #

const MIS_ =
    in(
        Set(Nucleus.GMT.read(joinpath(DA, "c2.all.v7.1.symbols.gmt"))["COLLER_MYC_TARGETS_UP"]),
    ).(MFE_)

# ---- #

const FE_X_SA_X_MSC = hcat(MSC_, MSC_ * 10.0, fill(0.8, lastindex(MFE_)))

# ---- #

const MSE_FE1_ = Nucleus.GMT.read(joinpath(DA, "h.all.v7.1.symbols.gmt"))

# ---- #

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
    ),
)

    @test isapprox(FeatureSetEnrichment._enrich!(al, MSC_, EX, MIS_, nothing), re; atol = 1e-12)

    # 43.375 μs (0 allocations: 0 bytes)
    # 37.583 μs (0 allocations: 0 bytes)
    # 164.917 μs (0 allocations: 0 bytes)
    # 186.417 μs (0 allocations: 0 bytes)
    # 326.042 μs (0 allocations: 0 bytes)
    # 325.916 μs (0 allocations: 0 bytes)
    #@btime FeatureSetEnrichment._enrich!($al, MSC_, EX, MIS_, nothing)

end

# ---- #

for al in AL_

    # 2.932 ms (108 allocations: 934.22 KiB)
    # 2.649 ms (108 allocations: 934.22 KiB)
    # 9.011 ms (108 allocations: 934.22 KiB)
    # 10.135 ms (108 allocations: 934.22 KiB)
    # 17.190 ms (108 allocations: 934.22 KiB)
    # 17.178 ms (108 allocations: 934.22 KiB)
    #@btime FeatureSetEnrichment.enrich($al, MFE_, MSC_, MFE1___)

end

# ---- #

for al in AL_

    # 9.282 ms (370 allocations: 5.51 MiB)
    # 8.426 ms (370 allocations: 5.51 MiB)
    # 27.542 ms (370 allocations: 5.51 MiB)
    # 30.772 ms (370 allocations: 5.51 MiB)
    # 52.249 ms (370 allocations: 5.51 MiB)
    # 52.131 ms (370 allocations: 5.51 MiB)
    #@btime FeatureSetEnrichment.enrich($al, MFE_, FE_X_SA_X_MSC, MFE1___)

end

# ---- #

const MSE_ = collect(keys(MSE_FE1_))

# ---- #

const MSA_ = ["Score", "Score x 10", "Constant"]

# ---- #

const AL = FeatureSetEnrichment.KS()

# ---- #

const SE_X_SA_X_EN = FeatureSetEnrichment.enrich(AL, MFE_, FE_X_SA_X_MSC, MFE1___)

# ---- #

FeatureSetEnrichment.plot(
    Nucleus.TE,
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
