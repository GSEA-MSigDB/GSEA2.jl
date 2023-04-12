using BioLab

using GSEA

# --------------------------------------------- #

TE = joinpath(tempdir(), "GSEA.test.small")

BioLab.Path.empty(TE)

# --------------------------------------------- #

DA = @__DIR__

sm = joinpath(DA, "small")

GSEA.metric_rank(
    joinpath(sm, "setting.json"),
    joinpath(DA, "target_x_sample_x_number.tsv"),
    joinpath(DA, "feature_x_sample_x_number.tsv"),
    joinpath(sm, "set_features.json"),
    TE,
)
