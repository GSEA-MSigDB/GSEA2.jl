function _tabulate_statistic(se_en, se_ra_, ou)

    se_ = collect(keys(se_en))

    en_ = collect(values(se_en))

    mkpath(ou)

    if isempty(se_ra_)

        gl_ = gla_ = fill(NaN, length(se_))

    else

        ra__ = [collect(values(se_ra)) for se_ra in se_ra_]

        gl_, gla_ = BioLab.Significance.get_p_value_and_adjust(en_, vcat(ra__...))

        se_x_ra_x_en = DataFrame("Set" => se_)

        insertcols!(se_x_ra_x_en, (string(id) => ra_ for (id, ra_) in enumerate(ra__))...)

        BioLab.Table.write(joinpath(ou, "set_x_random_x_enrichment.tsv"), se_x_ra_x_en)

    end

    se_x_st_x_nu = sort(
        DataFrame(
            "Set" => se_,
            "Enrichment" => en_,
            "Global p value" => gl_,
            "Adjusted global p value" => gla_,
        ),
        "Enrichment",
    )

    BioLab.Table.write(joinpath(ou, "set_x_statistic_x_number.tsv"), se_x_st_x_nu)

    se_x_st_x_nu

end
