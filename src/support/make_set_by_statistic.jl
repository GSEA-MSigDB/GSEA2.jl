function make_set_by_statistic(se_en, se_ra_, ou)

    if isempty(se_ra_)

        pv_ = ad_ = fill(NaN, length(se_en))

    else

        pv_ = []

        for (se, en) in se_en

            ra_ = []

            for se_ra in se_ra_

                ra = se_ra[se]

                if sign(en) == sign(ra)

                    push!(ra_, ra)

                end

            end

            if en < 0

                si = "<"

            else

                si = ">"

            end

            push!(pv_, get_p_value(en, ra_, si))

        end

        ad_ = adjust_p_value(pv_)

    end

    fl_se_st = sort(
        DataFrame(
            "Set" => collect(keys(se_en)),
            "Enrichment" => collect(values(se_en)),
            "P-value" => pv_,
            "Q-value" => ad_,
        ),
        "Enrichment",
    )

    mkpath(ou)

    table_write(joinpath(ou, "set_by_statistic.tsv"), fl_se_st)

    fl_se_st

end
