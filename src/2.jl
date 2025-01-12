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
