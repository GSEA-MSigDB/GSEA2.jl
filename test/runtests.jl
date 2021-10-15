using DataFrames

using GSEA
using Kwat

#
dw = joinpath(homedir(), "Downloads", "")

#
fe_, sc_, fe1_ = make_benchmark("card AK")

in_ = Kwat.vector.check_in(fe_, fe1_)

score_set(fe_, sc_, fe1_, in_)

score_set(fe_, sc_, fe1_; si = false, pa = joinpath(dw, "old_algorithm.png"))

score_set_new(fe_, sc_, fe1_; pa = joinpath(dw, "new_algorithm.jpeg"))

#
fe_, sc_, fe1_ = make_benchmark("myc")

in_ = Kwat.vector.check_in(fe_, fe1_)

score_set(fe_, sc_, fe1_, in_)

score_set(fe_, sc_, fe1_)

score_set_new(fe_, sc_, fe1_)

#
sc_fe_sa = DataFrame(
    "Feature" => fe_,
    "Score" => sc_,
    "Score x 10" => sc_ * 10,
    "Constant" => fill(0.8, length(fe_)),
)

se_fe_ = Kwat.gmt.read(joinpath("data", "h.all.v7.1.symbols.gmt"))

score_set(fe_, sc_, se_fe_)

score_set(sc_fe_sa, se_fe_; n_jo = 1)

#try_method(fe_, sc_, fe1_)
