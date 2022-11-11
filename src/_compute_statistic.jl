function _compute_statistic(se_en, se_ra__, di)

    #
    se_ = collect(keys(se_en))

    en_ = collect(values(se_en))

    #
    mkpath(di)

    #
    if isempty(se_ra__)

        #
        gl_ = gla_ = fill(NaN, length(se_))

    else

        #
        ra__ = [collect(values(se_ra_)) for se_ra_ in se_ra__]

        #
        gl_, gla_ = OnePiece.Significance.get_p_value_and_adjust(en_, vcat(ra__...))

        #
        se_x_ra_x_en = DataFrame("Set" => se_)

        insertcols!(se_x_ra_x_en, (id => ra_ for (id, ra_) in enumerate(ra__))...)

        OnePiece.Table.write(joinpath(di, "set_x_random_x_enrichment.tsv"), se_x_ra_x_en)

    end

    #
    se_x_st_x_nu = sort(
        DataFrame(
            "Set" => se_,
            "Enrichment" => en_,
            "Global p value" => gl_,
            "Adjusted global p value" => gla_,
        ),
        "Enrichment",
    )

    OnePiece.Table.write(joinpath(di, "set_x_statistic_x_number.tsv"), se_x_st_x_nu)

    se_x_st_x_nu

end
