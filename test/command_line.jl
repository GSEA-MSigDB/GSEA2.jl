const JS = joinpath(DA, "set_features.json")

const OU = joinpath(homedir(), "Downloads")

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
