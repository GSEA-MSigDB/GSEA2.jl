function plot_mountain(en_se_sa, n_ex, se_, sc_fe_sa, se_fe_, sy_ar, ou)

    pop!(sy_ar, :n_jo)

    for se in vcat(en_se_sa[1:n_ex, 1], se_)

        score_set(
            sc_fe_sa[!, 1],
            sc_fe_sa[!, 2],
            se_fe_[se];
            title_text = se,
            ou = clean,
            sy_ar...,
        )

    end

end
