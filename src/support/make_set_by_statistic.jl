function make_set_by_statistic(se_en, ra__, ou)

    se_ = collect(keys(se_en))

    en_ = collect(values(se_en))

    if isempty(ra__)

        pv_ = ad_ = fill(NaN, length(se_))

    else

        pv_, ad_ = get_p_value_and_adjust(en_, vcat(ra__...))

    end

    fl_se_st = sort(
        DataFrame("Set" => se_, "Enrichment" => en_, "P-Value" => pv_, "Q-Value" => ad_),
        "Enrichment",
    )

    mkpath(ou)

    table_write(joinpath(ou, "set_by_statistic.tsv"), fl_se_st)

    fl_se_st

end
