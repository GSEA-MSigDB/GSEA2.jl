@inline function _absolute_exponentiate(sc, ex)

    ab = abs(sc)

    if !isone(ex)

        ab ^= ex

    end

    ab

end

@inline function _get_1_normalizer(sc_, ex, is_)

    n = lastindex(sc_)

    su1 = 0.0

    for id in 1:n

        if is_[id]

            su1 += _absolute_exponentiate(sc_[id], ex)

        end

    end

    n, 1 / su1

end

@inline function _get_0_1_normalizer(sc_, ex, is_)

    n = lastindex(sc_)

    su0 = su1 = 0.0

    for id in 1:n

        if is_[id]

            su1 += _absolute_exponentiate(sc_[id], ex)

        else

            su0 += 1.0

        end

    end

    n, 1 / su0, 1 / su1

end

@inline function _get_all_1_normalizer(sc_, ex, is_)

    n = lastindex(sc_)

    su = su1 = 0.0

    for id in 1:n

        ab = _absolute_exponentiate(sc_[id], ex)

        su += ab

        if is_[id]

            su1 += ab

        end

    end

    n, 1 / su, 1 / su1

end

@inline function _get_0_normalizer(noa, no1)

    1 / (1 / noa - 1 / no1)

end

@inline function _clip(nu)

    ep = eps()

    nu < ep ? ep : nu

end

function enrich(al, fe_, sc_::AbstractVector, fe1___; mi = 1, ex = 1)

    en_ = Vector{Float64}(undef, lastindex(fe1___))

    fe_id = Dict(fe => id for (id, fe) in enumerate(fe_))

    for (id, fe1_) in enumerate(fe1___)

        is_ = Omics.Dict.is_in(fe_id, fe1_)

        en_[id] = sum(is_) < mi ? NaN : _enrich!(al, sc_, ex, is_, nothing)

    end

    en_

end

function enrich(al, fe_, fe_x_sa_x_sc, fe1___; mi = 1, ex = 1)

    se_x_sa_x_en = Matrix{Float64}(undef, lastindex(fe1___), size(fe_x_sa_x_sc, 2))

    no_ = BitVector(undef, lastindex(fe_))

    @showprogress for (id, sc_) in enumerate(eachcol(fe_x_sa_x_sc))

        no_ .= .!isnan.(sc_)

        scn_ = view(sc_, no_)

        so_ = sortperm(scn_; rev = true)

        se_x_sa_x_en[:, id] = enrich(al, view(fe_, no_)[so_], scn_[so_], fe1___; mi, ex)

    end

    se_x_sa_x_en

end
