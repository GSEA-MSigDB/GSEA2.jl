using DataFrames

using Test

using BioLab

using GSEA

# --------------------------------------------- #

di = "outer_loop"

pip = mkpath(joinpath(di, "python", "input"))

ip = mkpath(joinpath(di, "input"))

ou = mkpath(joinpath(di, "output"))

BioLab.Path.empty(ip)

# BioLab.Path.empty(ou)

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

# GSEA.metric_rank(joinpath(di, "setting.json"), tst, tsf, jss, ou; feature2_x_index_x_random)

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

_sen, se_, _st_, se_x_st_x_nu = BioLab.DataFrame.separate(daj)

n = length(se_)

set_x_index_x_random = BioLab.Table.read(joinpath(ou, "set_x_index_x_random.tsv"))

_sen, _se_, _id_, se_x_id_x_ra =
    BioLab.DataFrame.separate(set_x_index_x_random[indexin(se_, set_x_index_x_random[!, 1]), :])

# --------------------------------------------- #

men_ = Vector{Float64}(undef, n)

mep_ = Vector{Float64}(undef, n)

for id in 1:n

    ne_ = Vector{Float64}()

    po_ = Vector{Float64}()

    for ra in se_x_id_x_ra[id, :]

        if ra < 0.0

            push!(ne_, ra)

        elseif 0.0 < ra

            push!(po_, ra)

        end

    end

    men_[id] = mean(ne_)

    mep_[id] = mean(po_)

end

# --------------------------------------------- #

en_ = se_x_st_x_nu[:, 1]

@test en_ == dap[!, 5]

# --------------------------------------------- #

idp = findfirst(0.0 < en for en in en_)

idn_ = 1:(idp - 1)

idp_ = idp:n

enn_ = Vector{Float64}(undef, n)

for id in idn_

    enn_[id] = -en_[id] / men_[id]

end

for id in idp_

    enn_[id] = en_[id] / mep_[id]

end

@test all(isapprox(enn_[id], dap[id, 4]; atol = 10^-3) for id in 1:n)

# --------------------------------------------- #

ran_ = Vector{Float64}()

rap_ = Vector{Float64}()

for id2 in 1:size(se_x_id_x_ra, 2)

    for id1 in 1:size(se_x_id_x_ra, 1)

        ra = se_x_id_x_ra[id1, id2]

        if ra < 0.0

            push!(ran_, -ra / men_[id1])

        elseif 0.0 < ra

            push!(rap_, ra / mep_[id1])

        end

    end

end

pvl_, adl_ = BioLab.Significance.get_p_value_and_adjust(
    BioLab.Significance.get_p_value_for_less,
    enn_[idn_],
    ran_,
)

pvm_, adm_ = BioLab.Significance.get_p_value_and_adjust(
    BioLab.Significance.get_p_value_for_more,
    enn_[idp_],
    rap_,
)

pv_ = vcat(pvl_, pvm_)

ad_ = vcat(adl_, adm_)

@test all(isapprox(pv_[id], dap[id, 6]; atol = 10^-3) for id in 1:n)

@test all(isapprox(ad_[id], dap[id, 7]; atol = 10^-3) for id in 1:n)
