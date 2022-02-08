"""
Run standard GSEA

# Arguments

  - `js`:
  - `se`:
  - `tst`:
  - `tsd`:
  - `ou`: output directory
"""
@cast function run_standard_gsea(js, se, tst, tsd, ou)

    ke_ar = dict_read(js)

    nu_ta_sa = table_read(tst)

    sc_fe_sa = table_read(tsd)

    fe_ = string.(sc_fe_sa[:, 1])

    sc_ = compare_with_target(
        BitVector(nu_ta_sa[1, :]),
        Matrix(sc_fe_sa[:, 2:end]),
        "signal_to_noise_ratio",
    )

    se_fe_ = read_set(se, ke_ar)

    pe = pop!(ke_ar, "pe")

    if pe == "label"

        println("Permuting labels to compute significance")

        n_pe = pop!(ke_ar, "n_pe")

        ke_ar = symbolize_key(ke_ar)

        se_en = score_set(fe_, sc_, se_fe_; ke_ar...)

        # TODO: set random generator

        if 0 < n_pe

            sh_ = copy(sc_)

            _se_ra = []

            for it in 1:n_pe

                println("  ", it, "/", n_pe)

                push!(_se_ra, score_set(fe_, shuffle!(sh_), se_fe_; ke_ar...))

            end

            pv_, ad_ = get_p_value_and_adjust(se_en, se_ra_)

        else

            pv_ = ad_ = fill(NaN, length(se_en))

        end

        fl_se_st = make_set_by_statistic(se_en, pv_, ad_, ou)

        println("Plotting")

    elseif pe == "set"

        n_pe = pop!(ke_ar, "n_pe")

        run_pre_rank_gsea(ke_ar, se_fe_, fe_, sc_, n_pe, ou)

    else

        error(pe)

    end

end
