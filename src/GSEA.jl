module GSEA

using Printf: @sprintf

using Omics

struct KS end

struct KSa end

struct KLioM end

struct KLioP end

struct KLi end

struct KLi1 end

function _exponentiate(sc, ex)

    ab = abs(sc)

    isone(ex) ? ab : ab^ex

end

function _get_0_1(sc_, ex, is_)

    s0 = s1 = 0.0

    for id in eachindex(sc_)

        if is_[id]

            s1 += _exponentiate(sc_[id], ex)

        else

            s0 += 1.0

        end

    end

    -inv(s0), inv(s1)

end

function _get_all_1(sc_, ex, is_)

    su = s1 = 0.0

    for id in eachindex(sc_)

        so = _exponentiate(sc_[id], ex)

        su += so

        if is_[id]

            s1 += so

        end

    end

    inv(su), inv(s1)

end

function _get_0(no, n1)

    inv(inv(no) - inv(n1))

end

function _get_1(sc_, ex, is_)

    su = 0.0

    for id in eachindex(sc_)

        if is_[id]

            su += _exponentiate(sc_[id], ex)

        end

    end

    inv(su)

end

function _clip(fl)

    ep = eps()

    fl < ep ? ep : fl

end

function _enrich!(::KS, sc_, ex, is_, mo_)

    n0, n1 = _get_0_1(sc_, ex, is_)

    mo = ea = em = 0.0

    for id in eachindex(sc_)

        mo += is_[id] ? _exponentiate(sc_[id], ex) * n1 : n0

        if !isnothing(mo_)

            mo_[id] = mo

        end

        ab = abs(mo)

        if ea < ab

            ea = ab

            em = mo

        end

    end

    em

end

function _enrich!(::KSa, sc_, ex, is_, mo_)

    n0, n1 = _get_0_1(sc_, ex, is_)

    mo = ar = 0.0

    for id in eachindex(sc_)

        ar += mo += is_[id] ? _exponentiate(sc_[id], ex) * n1 : n0

        if !isnothing(mo_)

            mo_[id] = mo

        end

    end

    ar / lastindex(sc_)

end

function _enrich!(::KLioM, sc_, ex, is_, mo_)

    no, n1 = _get_all_1(sc_, ex, is_)

    n0 = _get_0(no, n1)

    ri = r0 = r1 = eps()

    pr = p0 = p1 = ar = 0.0

    le = l0 = l1 = 1.0

    for id in eachindex(sc_)

        so = _exponentiate(sc_[id], ex)

        de = so * no

        if is_[id]

            d0 = 0.0

            d1 = so * n1

        else

            d0 = so * n0

            d1 = 0.0

        end

        ar +=
            mo =
                Omics.Information.get_antisymmetric_kullback_leibler_divergence(
                    r1 += d1,
                    r0 += d0,
                    ri += de,
                ) - Omics.Information.get_antisymmetric_kullback_leibler_divergence(
                    _clip(l1 -= p1),
                    _clip(l0 -= p0),
                    _clip(le -= pr),
                )

        if !isnothing(mo_)

            mo_[id] = mo

        end

        pr = de

        p0 = d0

        p1 = d1

    end

    ar / lastindex(sc_)

end

function _enrich!(::KLioP, sc_, ex, is_, mo_)

    no, n1 = _get_all_1(sc_, ex, is_)

    n0 = _get_0(no, n1)

    ri = r0 = r1 = eps()

    pr = p0 = p1 = ar = 0.0

    le = l0 = l1 = 1.0

    for id in eachindex(sc_)

        so = _exponentiate(sc_[id], ex)

        de = so * no

        if is_[id]

            d0 = 0.0

            d1 = so * n1

        else

            d0 = so * n0

            d1 = 0.0

        end

        ar +=
            mo =
                Omics.Information.get_symmetric_kullback_leibler_divergence(
                    r1 += d1,
                    r0 += d0,
                    ri += de,
                ) - Omics.Information.get_symmetric_kullback_leibler_divergence(
                    _clip(l1 -= p1),
                    _clip(l0 -= p0),
                    _clip(le -= pr),
                )

        if !isnothing(mo_)

            mo_[id] = mo

        end

        pr = de

        p0 = d0

        p1 = d1

    end

    ar / lastindex(sc_)

end

function _enrich!(::KLi, sc_, ex, is_, mo_)

    no, n1 = _get_all_1(sc_, ex, is_)

    ri = r1 = eps()

    pr = p1 = ar = 0.0

    le = l1 = 1.0

    for id in eachindex(sc_)

        so = _exponentiate(sc_[id], ex)

        de = so * no

        d1 = is_[id] ? so * n1 : 0.0

        ar +=
            mo = Omics.Information.get_antisymmetric_kullback_leibler_divergence(
                r1 += d1,
                _clip(l1 -= p1),
                ri += de,
                _clip(le -= pr),
            )

        if !isnothing(mo_)

            mo_[id] = mo

        end

        pr = de

        p1 = d1

    end

    ar / lastindex(sc_)

end

function _enrich!(::KLi1, sc_, ex, is_, mo_)

    n1 = _get_1(sc_, ex, is_)

    ri = r1 = eps()

    de = inv(lastindex(sc_))

    p1 = ar = 0.0

    le = 1.0 + de

    l1 = 1.0

    for id in eachindex(sc_)

        d1 = is_[id] ? _exponentiate(sc_[id], ex) * n1 : 0.0

        ar +=
            mo = Omics.Information.get_antisymmetric_kullback_leibler_divergence(
                r1 += d1,
                _clip(l1 -= p1),
                ri += de,
                _clip(le -= de),
            )

        if !isnothing(mo_)

            mo_[id] = mo

        end

        p1 = d1

    end

    ar / lastindex(sc_)

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

    is_ = map(in(Set(me_)), fe_)

    mo_ = Vector{Float64}(undef, uf)

    en = _enrich!(al, sc_, ex, is_, mo_)

    tr = Dict("mode" => "lines", "line" => Dict("width" => 0), "fill" => "tozeroy")

    ie_ = findall(<(0.0), sc_)

    ip_ = findall(>=(0.0), sc_)

    da_ = [
        merge(tr, Dict("y" => sc_[ie_], "x" => xc_[ie_], "fillcolor" => Omics.Color.BL)),
        merge(tr, Dict("y" => sc_[ip_], "x" => xc_[ip_], "fillcolor" => Omics.Color.RE)),
        Dict(
            "yaxis" => "y2",
            "y" => zeros(sum(is_)),
            "x" => xc_[is_],
            "mode" => "markers",
            "marker" => Dict(
                "symbol" => "line-ns",
                "size" => 24,
                "line" => Dict(
                    "width" => 1.32,
                    "color" => Omics.Color.hexify(Omics.Color.SG, 0.8),
                ),
            ),
        ),
        merge(tr, Dict("yaxis" => "y3", "y" => mo_, "x" => xc_, "fillcolor" => "#07fa07")),
    ]

    if typeof(al) == KS

        id_ = findall(in(extrema(mo_)), mo_)

        push!(
            da_,
            Dict(
                "yaxis" => "y3",
                "y" => mo_[id_],
                "x" => xc_[id_],
                "mode" => "markers",
                "marker" =>
                    Dict("size" => 32, "color" => Omics.Color.hexify(Omics.Color.YE, 0.96)),
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
                "yaxis" => Dict("domain" => (0, 0.24), "title" => Dict("text" => "$ns")),
                "yaxis2" => Dict(
                    "domain" => (0.25, 0.31),
                    "title" => Dict("text" => "Set"),
                    "tickvals" => (),
                ),
                "yaxis3" =>
                    Dict("domain" => (0.32, 1), "title" => Dict("text" => "Δ Enrichment")),
                "xaxis" => Dict(
                    "title" => Dict("text" => "$nf ($uf)"),
                    "showspikes" => true,
                    "spikesnap" => "cursor",
                    "spikemode" => "across",
                    "spikedash" => "solid",
                    "spikethickness" => 1,
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

function _is_in!(is_, an_id, a1_)
    @assert !any(is_)

    for a1 in a1_

        id = get(an_id, a1, nothing)

        if !isnothing(id)

            is_[id] = true

        end

    end

end

function enrich(al, fe_, sc_::AbstractVector, me___; mi = 1, ex = 1)

    en_ = Vector{Float64}(undef, lastindex(me___))

    is_ = falses(lastindex(fe_))

    fe_id = Dict(fe_[id] => id for id in eachindex(fe_))

    for id in eachindex(me___)

        _is_in!(is_, fe_id, me___[id])

        en_[id] = sum(is_) < mi ? NaN : _enrich!(al, sc_, ex, is_, nothing)

        is_[is_] .= false

    end

    en_

end

function enrich(al, fe_, sc, me___; mi = 1, ex = 1)

    us = size(sc, 2)

    en = Matrix{Float64}(undef, lastindex(me___), us)

    no_ = BitVector(undef, lastindex(fe_))

    for id in 1:us

        sc_ = sc[:, id]

        map!(sc -> !isnan(sc), no_, sc_)

        so_ = sc_[no_]

        id_ = sortperm(so_; rev = true)

        en[:, id] = enrich(al, fe_[no_][id_], so_[id_], me___; mi, ex)

    end

    en

end

#include("command_line.jl")

end
