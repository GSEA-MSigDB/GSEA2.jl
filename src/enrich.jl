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

struct KS end

struct KSa end

struct KLi1 end

struct KLi end

struct KLioM end

struct KLioP end

function _enrich!(::KS, sc_, ex, is_, mo_)

    n, no0, no1 = _get_0_1_normalizer(sc_, ex, is_)

    cu = et = eta = 0.0

    for id in 1:n

        if is_[id]

            cu += _absolute_exponentiate(sc_[id], ex) * no1

        else

            cu -= no0

        end

        if !isnothing(mo_)

            mo_[id] = cu

        end

        cua = abs(cu)

        if eta < cua

            et = cu

            eta = cua

        end

    end

    et

end

function _enrich!(::KSa, sc_, ex, is_, mo_)

    n, no0, no1 = _get_0_1_normalizer(sc_, ex, is_)

    cu = ar = 0.0

    for id in 1:n

        if is_[id]

            cu += _absolute_exponentiate(sc_[id], ex) * no1

        else

            cu -= no0

        end

        if !isnothing(mo_)

            mo_[id] = cu

        end

        ar += cu

    end

    ar / n

end

function _enrich!(::KLi1, sc_, ex, is_, mo_)

    n, no1 = _get_1_normalizer(sc_, ex, is_)

    ri = ri1 = eps()

    rid = 1 / n

    le = 1.0 + rid

    le1 = 1.0

    ar = pr1 = 0.0

    for id in 1:n

        ri1d = is_[id] ? _absolute_exponentiate(sc_[id], ex) * no1 : 0.0

        en = Omics.Information.get_antisymmetric_kullback_leibler_divergence(
            ri1 += ri1d,
            _clip(le1 -= pr1),
            ri += rid,
            _clip(le -= rid),
        )

        ar += en

        if !isnothing(mo_)

            mo_[id] = en

        end

        pr1 = ri1d

    end

    ar / n

end

function _enrich!(::KLi, sc_, ex, is_, mo_)

    n, noa, no1 = _get_all_1_normalizer(sc_, ex, is_)

    ri = ri1 = eps()

    le = le1 = 1.0

    ar = pr = pr1 = 0.0

    for id in 1:n

        ab = _absolute_exponentiate(sc_[id], ex)

        ri1d = is_[id] ? ab * no1 : 0.0

        rid = ab * noa

        en = Omics.Information.get_antisymmetric_kullback_leibler_divergence(
            ri1 += ri1d,
            _clip(le1 -= pr1),
            ri += rid,
            _clip(le -= pr),
        )

        ar += en

        if !isnothing(mo_)

            mo_[id] = en

        end

        pr = rid

        pr1 = ri1d

    end

    ar / n

end

function _enrich!(::KLioM, sc_, ex, is_, mo_)

    n, noa, no1 = _get_all_1_normalizer(sc_, ex, is_)

    no0 = _get_0_normalizer(noa, no1)

    ri = ri1 = ri0 = eps()

    le = le1 = le0 = 1.0

    ar = pr = pr1 = pr0 = 0.0

    for id in 1:n

        ab = _absolute_exponentiate(sc_[id], ex)

        rid = ab * noa

        ri += rid

        le -= pr

        if is_[id]

            ri1d = ab * no1

            ri0d = 0.0

        else

            ri1d = 0.0

            ri0d = ab * no0

        end

        ri1 += ri1d

        le1 -= pr1

        ri0 += ri0d

        le0 -= pr0

        en =
            Omics.Information.get_antisymmetric_kullback_leibler_divergence(ri1, ri0, ri) -
            Omics.Information.get_antisymmetric_kullback_leibler_divergence(
                _clip(le1),
                _clip(le0),
                _clip(le),
            )

        ar += en

        if !isnothing(mo_)

            mo_[id] = en

        end

        pr = rid

        pr1 = ri1d

        pr0 = ri0d

    end

    ar / n

end

function _enrich!(::KLioP, sc_, ex, is_, mo_)

    n, noa, no1 = _get_all_1_normalizer(sc_, ex, is_)

    no0 = _get_0_normalizer(noa, no1)

    ri = ri1 = ri0 = eps()

    le = le1 = le0 = 1.0

    ar = pr = pr1 = pr0 = 0.0

    for id in 1:n

        ab = _absolute_exponentiate(sc_[id], ex)

        rid = ab * noa

        ri += rid

        le -= pr

        if is_[id]

            ri1d = ab * no1

            ri0d = 0.0

        else

            ri1d = 0.0

            ri0d = ab * no0

        end

        ri1 += ri1d

        le1 -= pr1

        ri0 += ri0d

        le0 -= pr0

        en =
            Omics.Information.get_symmetric_kullback_leibler_divergence(ri1, ri0, ri) -
            Omics.Information.get_symmetric_kullback_leibler_divergence(
                _clip(le1),
                _clip(le0),
                _clip(le),
            )

        ar += en

        if !isnothing(mo_)

            mo_[id] = en

        end

        pr = rid

        pr1 = ri1d

        pr0 = ri0d

    end

    ar / n

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
