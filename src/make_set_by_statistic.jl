function make_set_by_statistic(
    se_en::Dict{String, Float64},
    pv_::Vector{Float64},
    ad_::Vector{Float64},
    ou::String,
)::DataFrame

    fl_se_st = sort(
        DataFrame(
            "Set" => collect(keys(se_en)),
            "Enrichment" => collect(values(se_en)),
            "P-Value" => pv_,
            "Q-Value" => ad_,
        ),
        "Enrichment",
    )

    TableAccess.write(joinpath(ou, "set_by_statistic.tsv"), fl_se_st)

    return fl_se_st

end
