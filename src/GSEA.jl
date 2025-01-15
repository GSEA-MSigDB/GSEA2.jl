module GSEA

using Printf: @sprintf

using Omics

struct KS end

struct KSa end

struct KLioM end

struct KLioP end

struct KLi end

struct KLi1 end

function _is_in(fe_, me_)

    map(in(Set(me_)), fe_)

end

function _is_in!(bi_, fe_id, me_)

    for me in me_

        id = get(fe_id, me, nothing)

        if !isnothing(id)

            bi_[id] = true

        end

    end

end

function _exponentiate(sc, ex)

    ma = abs(sc)

    isone(ex) ? ma : ma^ex

end

function _get_0_1(sc_, ex, bo_)

    s0 = s1 = 0.0

    for id in eachindex(sc_)

        if bo_[id]

            s1 += _exponentiate(sc_[id], ex)

        else

            s0 += 1.0

        end

    end

    -inv(s0), inv(s1)

end

function _get_all_1(sc_, ex, bo_)

    sa = s1 = 0.0

    for id in eachindex(sc_)

        ma = _exponentiate(sc_[id], ex)

        sa += ma

        if bo_[id]

            s1 += ma

        end

    end

    inv(sa), inv(s1)

end

function _get_0(na, n1)

    inv(inv(na) - inv(n1))

end

function _get_1(sc_, ex, bo_)

    s1 = 0.0

    for id in eachindex(sc_)

        if bo_[id]

            s1 += _exponentiate(sc_[id], ex)

        end

    end

    inv(s1)

end

function _enrich!(::KS, sc_, ex, bo_, mo_)

    n0, n1 = _get_0_1(sc_, ex, bo_)

    mo = be = bs = 0.0

    for id in eachindex(sc_)

        mo += bo_[id] ? _exponentiate(sc_[id], ex) * n1 : n0

        if !isnothing(mo_)

            mo_[id] = mo

        end

        cu = abs(mo)

        if be < cu

            be = cu

            bs = mo

        end

    end

    bs

end

function _enrich!(::KSa, sc_, ex, bo_, mo_)

    n0, n1 = _get_0_1(sc_, ex, bo_)

    mo = ar = 0.0

    for id in eachindex(sc_)

        ar += mo += bo_[id] ? _exponentiate(sc_[id], ex) * n1 : n0

        if !isnothing(mo_)

            mo_[id] = mo

        end

    end

    ar / lastindex(sc_)

end

const ON = 1.0 + 1e-13

function _enrich!(::KLioM, sc_, ex, bo_, mo_)

    na, n1 = _get_all_1(sc_, ex, bo_)

    n0 = _get_0(na, n1)

    ra = r0 = r1 = eps()

    pa = p0 = p1 = ar = 0.0

    la = l0 = l1 = ON

    for id in eachindex(sc_)

        ma = _exponentiate(sc_[id], ex)

        da = ma * na

        if bo_[id]

            d0 = 0.0

            d1 = ma * n1

        else

            d0 = ma * n0

            d1 = 0.0

        end

        ra += da

        r0 += d0

        r1 += d1

        la -= pa

        l0 -= p0

        l1 -= p1

        ar +=
            mo =
                Omics.Information.get_antisymmetric_kullback_leibler_divergence(
                    r1,
                    r0,
                    ra,
                ) -
                Omics.Information.get_antisymmetric_kullback_leibler_divergence(l1, l0, la)

        if !isnothing(mo_)

            mo_[id] = mo

        end

        pa = da

        p0 = d0

        p1 = d1

    end

    ar / lastindex(sc_)

end

function _enrich!(::KLioP, sc_, ex, bo_, mo_)

    na, n1 = _get_all_1(sc_, ex, bo_)

    n0 = _get_0(na, n1)

    ra = r0 = r1 = eps()

    pa = p0 = p1 = ar = 0.0

    la = l0 = l1 = ON

    for id in eachindex(sc_)

        ma = _exponentiate(sc_[id], ex)

        da = ma * na

        if bo_[id]

            d0 = 0.0

            d1 = ma * n1

        else

            d0 = ma * n0

            d1 = 0.0

        end

        ra += da

        r0 += d0

        r1 += d1

        la -= pa

        l0 -= p0

        l1 -= p1

        ar +=
            mo =
                Omics.Information.get_symmetric_kullback_leibler_divergence(r1, r0, ra) -
                Omics.Information.get_symmetric_kullback_leibler_divergence(l1, l0, la)

        if !isnothing(mo_)

            mo_[id] = mo

        end

        pa = da

        p0 = d0

        p1 = d1

    end

    ar / lastindex(sc_)

end

function _enrich!(::KLi, sc_, ex, bo_, mo_)

    na, n1 = _get_all_1(sc_, ex, bo_)

    ra = r1 = eps()

    pa = p1 = ar = 0.0

    la = l1 = ON

    for id in eachindex(sc_)

        ma = _exponentiate(sc_[id], ex)

        da = ma * na

        d1 = bo_[id] ? ma * n1 : 0.0

        ra += da

        r1 += d1

        la -= pa

        l1 -= p1

        ar +=
            mo = Omics.Information.get_antisymmetric_kullback_leibler_divergence(
                r1,
                l1,
                ra,
                la,
            )

        if !isnothing(mo_)

            mo_[id] = mo

        end

        pa = da

        p1 = d1

    end

    ar / lastindex(sc_)

end

function _enrich!(::KLi1, sc_, ex, bo_, mo_)

    uf = lastindex(sc_)

    n1 = _get_1(sc_, ex, bo_)

    ra = r1 = eps()

    da = inv(uf)

    p1 = ar = 0.0

    la = ON + da

    l1 = ON

    for id in eachindex(sc_)

        d1 = bo_[id] ? _exponentiate(sc_[id], ex) * n1 : 0.0

        ra += da

        r1 += d1

        la -= da

        l1 -= p1

        ar +=
            mo = Omics.Information.get_antisymmetric_kullback_leibler_divergence(
                r1,
                l1,
                ra,
                la,
            )

        if !isnothing(mo_)

            mo_[id] = mo

        end

        p1 = d1

    end

    ar / uf

end

function _get_extreme(mo_)

    mi, ma = extrema(mo_)

    ai = abs(mi)

    aa = abs(ma)

    if ai == aa

        (mi, ma)

    elseif aa < ai

        (mi,)

    else

        (ma,)

    end

end

function plot(
    ht,
    al,
    fe_,
    sc_,
    me_;
    ex = 1.0,
    nf = "Feature",
    ns = "Score",
    nl = "Low",
    nh = "High",
    la = Dict{String, Any}(),
)

    uf = lastindex(fe_)

    xc_ = collect(1:uf)

    bo_ = _is_in(fe_, me_)

    mo_ = Vector{Float64}(undef, uf)

    en = _enrich!(al, sc_, ex, bo_, mo_)

    tr = Dict("mode" => "lines", "line" => Dict("width" => 0), "fill" => "tozeroy")

    ie_ = findall(<(0.0), sc_)

    ip_ = findall(>=(0.0), sc_)

    da_ = [
        merge(tr, Dict("y" => sc_[ie_], "x" => xc_[ie_], "fillcolor" => Omics.Color.BL)),
        merge(tr, Dict("y" => sc_[ip_], "x" => xc_[ip_], "fillcolor" => Omics.Color.RE)),
        Dict(
            "yaxis" => "y2",
            "y" => zeros(sum(bo_)),
            "x" => xc_[bo_],
            "mode" => "markers",
            "marker" => Dict(
                "symbol" => "line-ns",
                "size" => 24,
                "line" => Dict(
                    "width" => 2,
                    "color" => Omics.Color.hexify(Omics.Color.SG, 0.72),
                ),
            ),
        ),
        merge(tr, Dict("yaxis" => "y3", "y" => mo_, "x" => xc_, "fillcolor" => "#07fa07")),
    ]

    if typeof(al) == KS

        id_ = findall(in(_get_extreme(mo_)), mo_)

        push!(
            da_,
            Dict(
                "yaxis" => "y3",
                "y" => mo_[id_],
                "x" => xc_[id_],
                "mode" => "markers",
                "marker" =>
                    Dict("size" => 32, "color" => Omics.Color.hexify(Omics.Color.HU, 0.72)),
            ),
        )

    end

    an = Dict(
        "y" => 0,
        "font" => Dict("size" => 16),
        "borderpad" => 4.8,
        "borderwidth" => 2.64,
        "bordercolor" => Omics.Color.LI,
        "showarrow" => false,
    )

    ex = uf * 0.008

    Omics.Plot.plot(
        ht,
        da_,
        Omics.Dic.merg(
            Dict(
                "showlegend" => false,
                "yaxis" => Dict("domain" => (0, 0.24), "title" => Dict("text" => ns)),
                "yaxis2" => Dict(
                    "domain" => (0.248, 0.32),
                    "title" => Dict("text" => "Set"),
                    "tickvals" => (),
                ),
                "yaxis3" =>
                    Dict("domain" => (0.328, 1), "title" => Dict("text" => "Δ Enrichment")),
                "xaxis" => Dict(
                    "title" => Dict("text" => Omics.Strin.coun(uf, nf)),
                    "showspikes" => true,
                    "spikesnap" => "cursor",
                    "spikemode" => "across",
                    "spikedash" => "solid",
                    "spikethickness" => 0.8,
                    "spikecolor" => "#000000",
                ),
                "annotations" => (
                    Dict(
                        "yref" => "paper",
                        "xref" => "paper",
                        "y" => 1.064,
                        "text" => "Enrichment = <b>$(@sprintf "%.4g" en)</b>",
                        "font" => Dict("size" => 24, "color" => Omics.Color.BR),
                        "showarrow" => false,
                    ),
                    merge(
                        an,
                        Dict(
                            "x" => 1.0 - ex,
                            "xanchor" => "right",
                            "text" => nh,
                            "font" => Dict("color" => Omics.Color.RE),
                        ),
                    ),
                    merge(
                        an,
                        Dict(
                            "x" => uf + ex,
                            "xanchor" => "left",
                            "text" => nl,
                            "font" => Dict("color" => Omics.Color.BL),
                        ),
                    ),
                ),
            ),
            la,
        ),
    )

end

function enrich(al, fe_, sc_::AbstractVector, me___; um = 1, ex = 1)

    en_ = Vector{Float64}(undef, lastindex(me___))

    bi_ = falses(lastindex(fe_))

    fe_id = Dict(fe_[id] => id for id in eachindex(fe_))

    for id in eachindex(me___)

        _is_in!(bi_, fe_id, me___[id])

        en_[id] = sum(bi_) < um ? NaN : _enrich!(al, sc_, ex, bi_, nothing)

        bi_[bi_] .= false

    end

    en_

end

function enrich(al, fe_, sc, me___; um = 1, ex = 1)

    us = size(sc, 2)

    en = Matrix{Float64}(undef, lastindex(me___), us)

    bi_ = BitVector(undef, lastindex(fe_))

    for id in 1:us

        sc_ = sc[:, id]

        map!(!isnan, bi_, sc_)

        so_ = sc_[bi_]

        id_ = sortperm(so_; rev = true)

        en[:, id] = enrich(al, fe_[bi_][id_], so_[id_], me___; um, ex)

    end

    en

end

#include("command_line.jl")

end
