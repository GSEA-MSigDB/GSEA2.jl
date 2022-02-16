function make_set_x_statistic(se_en, se_ra_, ou)

    if isempty(se_ra_)

        no_ = pv_ = ad_ = fill(NaN, length(se_en))

    else

        no_ = []

        pv_ = []

        for (se, en) in se_en

            si = sign(en)

            ra_ = []

            for se_ra in se_ra_

                ra = se_ra[se]

                if sign(ra) == si

                    push!(ra_, ra)

                end

            end

            push!(no_, en / abs(mean(ra_)))

            push!(pv_, OnePiece.significance.get_p_value(en, ra_, si))

        end

        ad_ = OnePiece.significance.adjust_p_value(pv_)

    end

    fl_se_st = sort(
        DataFrame(
            "Set" => collect(keys(se_en)),
            "Enrichment" => collect(values(se_en)),
            "Normalized enrichment" => no_,
            "P-value" => pv_,
            "Q-value" => ad_,
        ),
        2,
    )

    mkpath(ou)

    OnePiece.table.write(joinpath(ou, "float.set_x_statistic.tsv"), fl_se_st)

    fl_se_st

end
