function plot_mountain(fl_se_st, n_ex, se_, fe_, sc_, se_fe_, sy_ar, ou)

    pop!(sy_ar, :n_jo)

    se2_ = fl_se_st[!, 1]

    di = joinpath(ou, "plot")

    mkpath(di)

    for se in vcat(se2_[1:n_ex], se2_[(end + 1 - n_ex):end], se_)

        score_set(
            fe_,
            sc_,
            se_fe_[se];
            title_text = se,
            ou = joinpath(di, string(clean(se, pr = false), ".html")),
            sy_ar...,
        )

    end

end
