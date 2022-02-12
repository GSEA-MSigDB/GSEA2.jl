OU = "set_by_statistic.tsv"

function make_set_by_statistic(se_en, pv_, ad_, ou)

    fl_se_st = sort(
        DataFrame(
            "Set" => collect(keys(se_en)),
            "Enrichment" => collect(values(se_en)),
            "P-Value" => pv_,
            "Q-Value" => ad_,
        ),
        "Enrichment",
    )

    mkpath(ou)

    table_write(joinpath(ou, OU), fl_se_st)

    fl_se_st

end
