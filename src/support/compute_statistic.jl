function compute_statistic(se_en, _se_ra, ou)

    se_ = collect(keys(se_en))

    en_ = collect(values(se_en))

    if isempty(_se_ra)

        gl_ = gla_ = fill(NaN, length(se_))

    else

        gl_, gla_ =
            OnePiece.significance.get_p_value_and_adjust(en_, vcat(collect.(values.(_se_ra))...))

    end

    se_x_st_x_nu = sort(
        DataFrame(
            "Set" => se_,
            "Enrichment" => en_,
            "Global pvalue" => gl_,
            "Adjusted global pvalue" => gla_,
        ),
        "Enrichment",
    )

    mkpath(ou)

    OnePiece.table.write(joinpath(ou, "set_x_statistic_x_number.tsv"), se_x_st_x_nu)

    se_x_st_x_nu

end
