using DataFrames

using Test

using BioLab

using GSEA

# --------------------------------------------- #

di = "test_outer_loop"

ip = mkpath(joinpath(di, "input"))

ou = mkpath(joinpath(di, "output"))

BioLab.Path.empty(ou)

# --------------------------------------------- #

fe_x_sa_x_sc = BioLab.GCT.read(
    "/Users/kwat/craft/GSEA_KS_run_all_files/input/Coller_et_al_gene_exp_preproc.gct",
)

sa_ = names(fe_x_sa_x_sc)[2:end]

rename!(fe_x_sa_x_sc, vcat("Gene", sa_))

fe = joinpath(ip, "feature_x_sample_x_score.tsv")

BioLab.Table.write(fe, fe_x_sa_x_sc)

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

ta = joinpath(ip, "target_x_sample_x_number.tsv")

BioLab.Table.write(ta, ta_x_sa_x_nu)

# --------------------------------------------- #

se_ge_ = BioLab.GMT.read((
    "/Users/kwat/craft/GSEA_KS_run_all_files/input/h.all.v2022.1.Hs.symbols.gmt",
    "/Users/kwat/craft/GSEA_KS_run_all_files/input/c2.all.v2022.1.Hs.symbols.gmt",
))

st = joinpath(ip, "set_genes.json")

BioLab.Dict.write(st, se_ge_)

# --------------------------------------------- #

GSEA.metric_rank(joinpath(ip, "setting.json"), ta, fe, st, ou)
