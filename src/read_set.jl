function read_set(
    se_fe_::String,
    ke_ar::Dict{String, Any},
)::Dict{String, Vector{String}}

    return select(gmt_read(se_fe_), pop!(ke_ar, "mi"), pop!(ke_ar, "ma"))

end

export read_set
