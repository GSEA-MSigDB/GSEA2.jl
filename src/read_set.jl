function read_set(
    gm::String,
    ke_ar::Dict{String, T} where {T <: Any},
)::Dict{String, Vector{String}}

    return select_set(GMTAccess.read(gm), pop!(ke_ar, "mi"), pop!(ke_ar, "ma"))

end
