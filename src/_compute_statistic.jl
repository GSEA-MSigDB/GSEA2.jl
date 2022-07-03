function _compute_statistic(se_en, se_ra__, ou)

    se_ = collect(keys(se_en))

    en_ = collect(values(se_en))

    mkpath(ou)

    if isempty(se_ra__)

        gl_ = gla_ = fill(NaN, length(se_))

    else

        ra__ = collect.(values.(se_ra__))

        gl_, gla_ = OnePiece.significance.get_p_value_and_adjust(en_, vcat(ra__...))

        se_x_ra_x_en = DataFrame("Set" => se_)

        insertcols!(se_x_ra_x_en, (string(id) => ra_ for (id, ra_) in enumerate(ra__))...)

        OnePiece.table.write(joinpath(ou, "set_x_random_x_enrichment.tsv"), se_x_ra_x_en)

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

    OnePiece.table.write(joinpath(ou, "set_x_statistic_x_number.tsv"), se_x_st_x_nu)

    se_x_st_x_nu

end
