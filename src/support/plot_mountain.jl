function plot_mountain(se_x_st_x_nu, n_ex, pl_, fe_, sc_, se_fe_, sy_ar, ou)

    se_ = se_x_st_x_nu[!, 1]

    di = joinpath(ou, "plot")

    mkpath(di)

    pop!(sy_ar, :n_jo)

    for se in vcat(se_[1:n_ex], se_[(end - n_ex + 1):end], pl_)

        OnePiece.feature_set_enrichment.score_set(
            fe_,
            sc_,
            se_fe_[se],
            title_text = se,
            ou = joinpath(di, "$(OnePiece.path.clean(se, pr = false)).html"),
            sy_ar...,
        )

    end

end
