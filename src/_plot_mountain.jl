function _plot_mountain(se_x_st_x_nu, fe, sc, n_ex, pl_, al, fe_, sc_, se_fe_, sy_ar, di)

    n_se = size(se_x_st_x_nu, 1)

    n_ex = min(n_ex, n_se)

    # TODO: Try `1:2`.
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

    pl = mkpath(joinpath(di, "plot"))

    pop!(sy_ar, :n_jo)

    for se in pl_

        BioLab.FeatureSetEnrichment.score_set(
            al,
            fe_,
            sc_,
            se_fe_[se];
            title_text = se,
            fe = fe,
            sc = sc,
            ht = joinpath(pl, "$(BioLab.Path.clean(se)).html"),
            sy_ar...,
        )

    end

end
