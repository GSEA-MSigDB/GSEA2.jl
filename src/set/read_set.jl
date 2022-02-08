function read_set(
    pa::String,
    ke_ar::Dict{String, T} where {T <: Any},
)::Dict{String, Vector{String}}

    ex = splitext(pa)[2]

    if ex == ".gmt"

        se_ge_ = GMTAccess.read(pa)

    elseif ex == ".json"

        se_ge_ = convert(Dict{String, Vector{String}}, DictExtension.read(pa))

    else

        error(ex)

    end

    return select_set(se_ge_, pop!(ke_ar, "mi"), pop!(ke_ar, "ma"))

end
