function _plot_mountain(se_x_st_x_nu, n_ex, pl_, al, fe_, sc_, se_fe_, sy_ar, ou)

    n_se = size(se_x_st_x_nu, 1)

    n_ex = minimum([n_ex, n_se])

    co_ = [1, 2]

    for ro in 1:n_ex

        se, en = se_x_st_x_nu[ro, co_]

        if en <= 0 && !(se in pl_)

            push!(pl_, se)

        end

    end

    for ro in n_se:-1:(n_se - n_ex + 1)

        se, en = se_x_st_x_nu[ro, co_]

        if 0 <= en && !(se in pl_)

            push!(pl_, se)

        end

    end

    di = joinpath(ou, "plot")

    mkpath(di)

    pop!(sy_ar, :n_jo)

    if al == "cidac"

        fu = OnePiece.feature_set_enrichment.score_set_new

    elseif al == "ks"

        fu = OnePiece.feature_set_enrichment.score_set

    end

    for se in pl_

        fu(
            fe_,
            sc_,
            se_fe_[se];
            title_text = se,
            fe = "Gene",
            ou = joinpath(di, "$(OnePiece.path.clean(se, pr = false)).html"),
            sy_ar...,
        )

    end

end
