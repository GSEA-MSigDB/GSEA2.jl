using DataFrames: DataFrame, names

using Kwat.math: get_center
using Kwat.vector: check_in, sort_like

function _sum_1_absolute_and_0_count(
    sc_::Vector{Float64},
    in_::Vector{Float64},
)::Tuple{Float64, Float64}

    su1 = 0.0

    su0 = 0.0

    for id in 1:length(sc_)

        if in_[id] == 1.0

            nu = sc_[id]

            if nu < 0.0

                nu = -nu

            end

            su1 += nu

        else

            su0 += 1.0

        end

    end

    return su1, su0

end

function score_set(
    fe_::Vector{String},
    sc_::Vector{Float64},
    fe1_::Vector{String},
    in_::Vector{Float64};
    we::Float64 = 1.0,
    al::String = "ks",
    pl::Bool = true,
    ke_ar...,
)::Float64

    n_fe = length(fe_)

    en = 0.0

    en_ = Vector{Float64}(undef, n_fe)

    ex = 0.0

    exa = 0.0

    ar = 0.0

    su1, su0 = _sum_1_absolute_and_0_count(sc_, in_)

    de = 1.0 / su0

    @inbounds @fastmath @simd for id in n_fe:-1:1

        if in_[id] == 1.0

            sc = sc_[id]

            if sc < 0.0

                sc = -sc

            end

            en += sc / su1

        else

            en -= de

        end

        if pl

            en_[id] = en

        end

        ar += en

        if en < 0.0

            ena = -en

        else

            ena = en

        end

        if exa < ena

            exa = ena

            ex = en

        end

    end

    if al == "ks"

        en = ex

    elseif al == "auc"

        en = ar / Float64(n_fe)

    end

    if pl

        plot_mountain(fe_, sc_, in_, en_, en; ke_ar...)

    end

    return en

end

function score_set(
    fe_::Vector{String},
    sc_::Vector{Float64},
    fe1_::Vector{String};
    we::Float64 = 1.0,
    al::String = "ks",
    pl::Bool = true,
    ke_ar...,
)::Float64

    return score_set(
        fe_,
        sc_,
        fe1_,
        check_in(fe_, fe1_);
        we = we,
        al = al,
        pl = pl,
        ke_ar...,
    )

end

function score_set(
    fe_::Vector{String},
    sc_::Vector{Float64},
    se_fe_::Dict{String, Vector{String}};
    we::Float64 = 1.0,
    al::String = "ks",
)::Dict{String, Float64}

    if length(se_fe_) < 10

        ch = fe_

    else

        ch = Dict(fe => id for (fe, id) in zip(fe_, 1:length(fe_)))

    end

    se_en = Dict{String, Float64}()

    for (se, fe1_) in se_fe_

        se_en[se] = score_set(
            fe_,
            sc_,
            fe1_,
            check_in(ch, fe1_);
            we = we,
            al = al,
            pl = false,
        )

    end

    return se_en

end

function score_set(
    sc_fe_sa::DataFrame,
    se_fe_::Dict{String, Vector{String}};
    we::Float64 = 1.0,
    al::String = "ks",
    n_jo::Int64 = 1,
)::DataFrame

    fe_ = sc_fe_sa[!, 1]

    en_se_sa = DataFrame("Set" => collect(keys(se_fe_)))

    for sa in names(sc_fe_sa)[2:end]

        go_ = findall(!ismissing, sc_fe_sa[!, sa])

        sc_, fe_ = sort_like([sc_fe_sa[go_, sa], fe_[go_]])

        if in(al, ["ks", "auc"])

            se_en = score_set(fe_, sc_, se_fe_; we = we, al = al)

        elseif al == "js"

            se_en = score_set_new(fe_, sc_, se_fe_)

        end

        en_se_sa[!, sa] = collect(se_en[se] for se in en_se_sa[!, "Set"])

    end

    return en_se_sa

end

export score_set
