function make_set_by_statistic(se_en, ra__, ou)

    se_ = collect(keys(se_en))

    en_ = collect(values(se_en))

    if isempty(ra__)

        pv_ = ad_ = fill(NaN, length(se_))

    else

        ra_ = vcat(ra__...)

        enn_ = en_ .< 0

        ran_ = ra_ .< 0

        enp_ = .!enn_

        rap_ = .!ran_

        pvn_, adn_ = get_p_value_and_adjust(en_[enn_], ra_[ran_], "<")

        pvp_, adp_ = get_p_value_and_adjust(en_[enp_], ra_[rap_], ">")

        se_ = vcat(se_[enn_], se_[enp_])

        en_ = vcat(en_[enn_], en_[enp_])

        pv_ = vcat(pvn_, pvp_)

        ad_ = vcat(adn_, adp_)

    end

    fl_se_st = sort(
        DataFrame("Set" => se_, "Enrichment" => en_, "P-value" => pv_, "Q-value" => ad_),
        "Enrichment",
    )

    mkpath(ou)

    table_write(joinpath(ou, "set_by_statistic.tsv"), fl_se_st)

    fl_se_st

end
