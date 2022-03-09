function make_set_x_statistic(se_en, _se_ra, ou)

    se_ = collect(keys(se_en))

    en_ = collect(values(se_en))

    if isempty(_se_ra)

        no_ = lo_ = loa_ = gl_ = gla_ = fill(NaN, length(se_en))

    else

        no_ = []

        lo_ = []

        for (se, en) in se_en

            si = sign(en)

            ra_ = []

            for se_ra in _se_ra

                ra = se_ra[se]

                if sign(ra) == si

                    push!(ra_, ra)

                end

            end

            push!(no_, en / abs(mean(ra_)))

            push!(lo_, OnePiece.significance.get_p_value(en, ra_, si))

        end

        loa_ = OnePiece.significance.adjust_p_value(lo_)

        gl_, gla_ =
            OnePiece.significance.get_p_value_and_adjust(en_, vcat(collect.(values.(_se_ra))...))

    end

    se_x_st = sort(
        DataFrame(
            "Set" => se_,
            "Enrichment" => en_,
            "Gene-set-size-normalized enrichment" => no_,
            "Local pvalue" => lo_,
            "Adjusted local pvalue" => loa_,
            "Global pvalue" => gl_,
            "Adjusted global pvalue" => gla_,
        ),
        "Enrichment",
    )

    mkpath(ou)

    OnePiece.table.write(joinpath(ou, "set_x_statistic_x_number.tsv"), se_x_st)

    se_x_st

end
