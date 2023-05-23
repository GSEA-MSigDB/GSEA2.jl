using BioLab

using GSEA

# ---- #

TE = joinpath(tempdir(), "GSEA.test.small")

BioLab.Path.empty(TE)

# ---- #

DA = @__DIR__

di = joinpath(DA, "small")

GSEA.metric_rank(
    joinpath(di, "setting.json"),
    joinpath(DA, "target_x_sample_x_number.tsv"),
    joinpath(DA, "feature_x_sample_x_number.tsv"),
    joinpath(di, "set_features.json"),
    TE,
)
