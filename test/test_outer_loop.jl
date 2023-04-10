using DataFrames

using Test

using BioLab

using GSEA

# --------------------------------------------- #

di = "test_outer_loop"

pip = mkpath(joinpath(di, "python", "input"))

ip = mkpath(joinpath(di, "input"))

ou = mkpath(joinpath(di, "output"))

BioLab.Path.empty(ip)

BioLab.Path.empty(ou)

# --------------------------------------------- #

_fen, fe_, sa_, fe_x_sa_x_nu =
    BioLab.DataFrame.separate(BioLab.GCT.read(joinpath(pip, "Coller_et_al_gene_exp_preproc.gct")))

BioLab.Matrix.apply_by_row!(BioLab.Normalization.normalize_with_0!, fe_x_sa_x_nu)

clamp!(fe_x_sa_x_nu, -3, 3)

tsf = joinpath(ip, "feature_x_sample_x_number.tsv")

BioLab.Table.write(tsf, BioLab.DataFrame.make("Feature", fe_, sa_, fe_x_sa_x_nu))

# --------------------------------------------- #

ta_ = replace(
    split(readlines(joinpath(pip, "Coller_et_al_phen.cls"))[3], ' '; keepempty = false),
    "cntrl" => 0,
    "myc" => 1,
)

target_x_sample_x_number = DataFrame(permutedims(vcat("Control vs MYC", ta_)), vcat("Target", sa_))

tst = joinpath(ip, "target_x_sample_x_number.tsv")

BioLab.Table.write(tst, target_x_sample_x_number)

# --------------------------------------------- #

se_fe_ = BioLab.GMT.read((
    joinpath(pip, "h.all.v2022.1.Hs.symbols.gmt"),
    joinpath(pip, "c2.all.v2022.1.Hs.symbols.gmt"),
))

jss = joinpath(ip, "set_features.json")

BioLab.Dict.write(jss, se_fe_)

# --------------------------------------------- #

feature2_x_index_x_random = BioLab.Table.read(
    joinpath(di, "python", "txt", "KS_SUP example mean scaling_rand_perm_gene_scores.txt"),
)

GSEA.metric_rank(joinpath(di, "setting.json"), tst, tsf, jss, ou; feature2_x_index_x_random)

# --------------------------------------------- #

dap = BioLab.Table.read(
    joinpath(di, "python", "txt", "KS_SUP example mean scaling_gene_selection_scores.txt"),
)

daj = BioLab.Table.read(joinpath(ou, "feature_x_metric_x_score.tsv"))

idp = 2

idj = 2

digits = 3

transform!(dap, idp => co -> [round(nu; digits) for nu in co]; renamecols = false)

transform!(daj, idj => co -> [round(nu; digits) for nu in co]; renamecols = false)

sort!(dap, [idp, 1])

sort!(daj, [idj, 1])

@test isequal(dap[!, 1], daj[!, 1])

@test isequal(dap[!, idp], daj[!, idj])

# --------------------------------------------- #

dap = BioLab.Table.read(
    joinpath(di, "python", "txt", "KS_SUP example mean scaling_GSEA_results_table.txt"),
)

daj = BioLab.Table.read(joinpath(ou, "set_x_statistic_x_number.tsv"))

idp = 5

idj = 2

digits = 3

transform!(dap, idp => co -> [round(nu; digits) for nu in co]; renamecols = false)

transform!(daj, idj => co -> [round(nu; digits) for nu in co]; renamecols = false)

sort!(dap, [idp, 1])

sort!(daj, [idj, 1])

@test isequal(dap[!, 1], daj[!, 1])

@test isequal(dap[!, idp], daj[!, idj])

# --------------------------------------------- #

fe_x_id_x_ra
