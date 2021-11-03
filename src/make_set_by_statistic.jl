function make_set_by_statistic(
    se_en::Dict{String, Float64},
    pv_::Vector{Float64},
    ad_::Vector{Float64},
)::DataFrame

    return sort(
        DataFrame(
            "Set" => collect(keys(se_en)),
            "Enrichment" => collect(values(se_en)),
            "P-Value" => pv_,
            "Q-Value" => ad_,
        ),
        "Enrichment",
    )

end

export make_set_by_statistic
