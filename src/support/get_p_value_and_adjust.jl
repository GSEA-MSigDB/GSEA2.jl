function get_p_value_and_adjust(se_en, se_ra_)

    pv_ = []

    for (se, en) in se_en

        if en < 0

            si = "<"

        else

            si = ">"

        end

        push!(pv_, get_p_value(en, [se_ra[se] for se_ra in se_ra_], si))

    end

    pv_, adjust_p_value(pv_)

end
