function read_set(pa, ke_ar)

    ex = splitext(pa)[2]

    if ex == ".gmt"

        se_ge_ = gmt_read(pa)

    elseif ex == ".json"

        se_ge_ = dict_read(pa)

    else

        error(ex)

    end

    select_set(se_ge_, pop!(ke_ar, "mi"), pop!(ke_ar, "ma"))

end
