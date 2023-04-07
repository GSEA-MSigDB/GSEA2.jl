using DataFrames

using Test

using BioLab

using GSEA

# --------------------------------------------- #

di = mkpath("test_outer_loop")

# --------------------------------------------- #

fe_x_sa_x_sc = BioLab.GCT.read(
    "/Users/kwat/craft/GSEA_KS_run_all_files/input/Coller_et_al_gene_exp_preproc.gct",
)

sa_ = names(fe_x_sa_x_sc)[2:end]

rename!(fe_x_sa_x_sc, vcat("Gene", sa_))

BioLab.Table.write(joinpath(di, "feature_x_sample_x_score.tsv"), fe_x_sa_x_sc)

# --------------------------------------------- #

ta_ = replace(
    split(
        readlines("/Users/kwat/craft/GSEA_KS_run_all_files/input/Coller_et_al_phen.cls")[3],
        ' ';
        keepempty = false,
    ),
    "cntrl" => 0,
    "myc" => 1,
)

ta_x_sa_x_nu = DataFrame(permutedims(vcat("Control vs MYC", ta_)), vcat("Target", sa_))

BioLab.Table.write(joinpath(di, "target_x_sample_x_number.tsv"), ta_x_sa_x_nu)

# --------------------------------------------- #

se_ge_ =
    BioLab.GMT.read("/Users/kwat/craft/GSEA_KS_run_all_files/input/h.all.v2022.1.Hs.symbols.gmt")

BioLab.Dict.write(joinpath(di, "set_genes.json"), se_ge_)
