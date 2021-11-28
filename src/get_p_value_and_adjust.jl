function get_p_value_and_adjust(
    se_en::Dict{String, Float64},
    se_ra_::Vector{Dict{String, Float64}},
)::Tuple{Vector{Float64}, Vector{Float64}}

    pv_ = Vector{Float64}()

    for (se, en) in se_en

        if en < 0

            si = "<"

        else

            si = ">"

        end

        ra_ = [se_ra[se] for se_ra in se_ra_]

        push!(pv_, Significance.get_p_value(en, ra_, si))

    end

    return pv_, Significance.adjust_p_value(pv_)

end
